"""
  pybsseq
"""
import numpy as np
import pandas as pd
import argparse
import logging
import sys
import os, os.path
import scipy.stats as st
import allel
import subprocess
from methylpy.call_mc import *
from methylpy.DMRfind import DMRfind
import glob
import h5py
import re
import subprocess, shlex

log = logging.getLogger(__name__)

def die(msg):
  sys.stderr.write('Error: ' + msg + '\n')
  sys.exit(1)

def readVcf(inFile):
    bvcf = allel.read_vcf(inFile, samples = [0], fields = '*')
    ## Changed this from vcfnp to allel, check the rest of scripts
    #bvcf = vcfnp.variants(inFile, cache=True).view(np.recarray)
    #bvcfD = vcfnp.calldata_2d(inFile, cache=True).view(np.recarray)
    return(bvcf)

def BinomTest(per_n, n, p, alternative="greater"):
    if per_n < 1 and per_n > 0:
        tpVal = st.binom_test(per_n * n, n, p, alternative = alternative)
    else:
        tpVal = st.binom_test(per_n, n, p, alternative = alternative)
    return tpVal

def getConvRate(bsCHROM, bsCONTEXT, bsMethPer, bstC, chrs = "ChrC"):
  chrInd = np.where(bsCHROM == chrs)[0]
  contInd = chrInd[np.where((bsCONTEXT[chrInd] == 'CG') | (bsCONTEXT[chrInd] == 'CHG') | (bsCONTEXT[chrInd] == 'CHH'))[0]]
  chrMethPer = bsMethPer[contInd]
  chrDepth = bstC[contInd]
  conv_rate = np.nansum(np.multiply(chrMethPer, chrDepth))/np.nansum(chrDepth)
  return conv_rate

def callMPs(bsMethPer, bstC, error_rate, alternative="greater", window=300000):
  bsPval = np.zeros(0,dtype=float)
  npBinomTest = np.vectorize(BinomTest)
  log.info("running binomial test for indiviadual sites!")
  for i in range(0, len(bsMethPer), window):
    pVal = npBinomTest(bsMethPer[i:i+window], bstC[i:i+window], error_rate, alternative=alternative)
    bsPval = np.append(bsPval, pVal)
    log.info("progress: %s positions" % (i + window))
  return bsPval

def getconv_rate_allc(allcFile):
    ## Give allcfile for the Chromosome of unMethylatedControl
    bsbed = pd.read_table(allcFile)
    total_unconverted_c = np.sum(bsbed['mc_count'])
    total_bases = np.sum(bsbed['total'])
    non_conversion = total_unconverted_c / float(total_bases)
    conv_rate = 1 - non_conversion
    log.info("Conversion rate estimated from %s is: %s" % (allcFile, conv_rate))
    return conv_rate

def getConvRate_mpipleup(path_to_allc, sample_id, control="ChrC"):
    ## Copied this code from methylpy to get the Conversion rate and error rate
    #figure out non-conversion rate if it isn't given
    non_conversion = 0
    quality_base = 33
    mpileupFile = glob.glob(path_to_allc + sample_id + "*" + control + "*tsv")
    try:
        f = open(mpileupFile[0] ,'r')
    except IOError:
        die("mpileup file not present")
    unconverted_c = 0
    converted_c = 0
    total_phred = 0
    total_bases = 0
    total_unconverted_c =0
    total_converted_c = 0
    for line in f:
        line = line.rstrip()
        fields = line.split("\t")
        if fields[3] != "0":
            total_phred += sum([ord(i)-quality_base for i in fields[5]])
            total_bases += len(fields[5])
            if fields[2] == "C":
                unconverted_c = fields[4].count(".")
                converted_c = fields[4].count("T")
            elif fields[2] == "G":
                unconverted_c = fields[4].count(",")
                converted_c = fields[4].count("a")
            else:
                continue
            total_unconverted_c += unconverted_c
            total_converted_c += converted_c
    f.close()
    #compute pvalues for control genome. Have to read through it completely twice to do this
    min_pvalue = 1
    avg_qual = total_phred / total_bases
    seq_error = 10 ** (avg_qual / -10)
    if seq_error > 1.0 or seq_error < 0.0:
        die("One of your quality values corresponds to a sequence error rate of "+str(seq_error)+". These error rates have to be between 0 and 1. Are you sure you chose the correct CASAVA version?")
    non_conversion = total_unconverted_c / float(total_converted_c + total_unconverted_c)
    log.info("\tThe non-conversion rate is "+str(non_conversion*100)+"%")
    log.info("\tThe estimated sequencing error rate is: "+str(seq_error))
    conv_rate = 1 - non_conversion
    return (conv_rate, seq_error)

def writeBED(bsCHROM, bsPOS, bsCONTEXT, bstC, bsMethPer, bsPval, bsSTRAND, outBED):
  out = open(outBED, 'w')
  for i in len(bsPOS):
    out.write("%s\t%s\t%s\t%s:%s:%s\t%s\t%s\n" % (bsCHROM[i], bsPOS[i], bsPOS[i]+1, bsCONTEXT[i], bstC[i], bsMethPer[i], bsPval[i], bsSTRAND[i]))
  out.close()

def getMPsfromVCF(args):
  (bVCF, bvcfD) = readVCF(args['inVCF'])
  error_rate = getConvRate(bvcf.CHROM, bvcf.CX, bvcfD.BT[:,0], bvcfD.CV[:,0], chrs = "ChrC")
  log.info("Conversion rate: %s", 100 - error_rate * 100)
  ChrsNums = np.array(("Chr1","Chr2","Chr3","Chr4","Chr5"))
  MethInd = np.where((np.in1d(bvcf.CHROM, ChrsNums)) & (np.in1d(bvcf.CX, MethContext)))[0]
  log.info("Number of positions: %s", len(MethInd))
  bsPval = callMPs(bvcfD.BT[MethInd,0], bvcfD.CV[MethInd,0], error_rate, window=args['window'])
  bsSTRAND = np.core.defchararray.replace(np.core.defchararray.replace(bvcf.REF[MethInd], "C", "+"), "G", "-")
  log.info("writing MPs info in out bed file")
  writeBED(bvcf.CHROM[MethInd], bvcf.POS[MethInd], bvcf.CX[MethInd], bvcfD.CV[MethInd,0], bvcfD.BT[MethInd,0], bsPval, bsSTRAND, outBED)
  log.info("finished!")

def methylpy_callmcs(args):
    files = [args['inFile']]
    libraries = [1]
    sample = args['sample_id']
    f_reference = args['ref_fol'] + '_f'  ## for forward ref
    r_reference = args['ref_fol'] + '_r'
    reference_fasta = args['ref_fasta']
    uMeth = args['unMeth']
    mem = args['memory']
    procs = args['nt']
    run_methylation_pipeline(files,libraries,sample,f_reference,r_reference,reference_fasta,uMeth,sort_mem=mem,num_procs=procs,sig_cutoff=0.01,min_cov=3)

def getChrs_allc(allcFile):
    ## path to one allc file
    bsbed = pd.read_table(allcFile, nrows = 5)  ## just read top 100 lines
    inds = np.core.defchararray.rfind(np.array(bsbed['chr'], dtype="string"),'Chr')
    if np.unique(inds)[0] >= 0:
        chrs = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
        return chrs
    inds = np.core.defchararray.rfind(np.array(bsbed['chr'], dtype="string"),'chr')
    if np.unique(inds)[0] >= 0:
        chrs = ["chr1", "chr2", "chr3", "chr4", "chr5"]
        return chrs
    return ["1", "2", "3", "4", "5"]

def methylpy_dmrfind(args):
    region_dict = {}
    checkAllcFile = glob.glob(args['path_to_allc']+"/*_1.tsv")
    chrs = getChrs_allc(checkAllcFile[0])
    for chrom in chrs:
        region_dict[chrom] = [0, 5000000000]
    mc_type = args['mc_type'].split(",")
    samples = args['sample_id'].split(",")
    if args['sample_cat'] == '0':
        sample_cat=range(0, len(samples))
    else:   ## if samples given some categories
        sample_cat = args['sample_cat'].split(",")
    ### Check which of the samples given have the all c files
    dbacc=np.zeros(0)
    for eacc in checkAllcFile:
        accID = os.path.basename(eacc).replace("allc_", "").replace("_1.tsv", "")
        dbacc = np.append(dbacc, accID)
    reqAcc=np.array(samples)[np.where(np.in1d(np.array(samples), dbacc))[0]]
    log.info("Number of accessions present in allc: %s", len(reqAcc))
    log.info(reqAcc.tolist())
    DMRfind(mc_type=mc_type, region_dict=region_dict, samples=reqAcc.tolist(), path_to_allc=args['path_to_allc'], num_procs=args['nt'], save_result=args['outDMR'], min_cov=1)

def parseMC_class(mc_class):
    if mc_class[1].upper() == 'G':
        context = ["CG",0]
    elif mc_class[2].upper() == 'G':
        context = ["CHG",1]
    elif mc_class:
        context = ["CHH",2]
    return context

def outH5File(allcBed, chrpositions, outFile):
    h5file = h5py.File(outFile, 'w')
    num_lines = len(chrpositions)
    h5file.create_dataset('chrpositions', data=chrpositions, shape=(num_lines,),dtype='i4')
    ## Going through all the columns
    for i in range(allcBed.shape[1]):
        try:
            h5file.create_dataset(allcBed.columns[i], compression="gzip", data=np.array(allcBed[allcBed.columns[i]]), shape=(allcBed.shape[0],))
        except TypeError:
            h5file.create_dataset(allcBed.columns[i], compression="gzip", data=np.array(allcBed[allcBed.columns[i]]).tolist(), shape=(allcBed.shape[0],))
    h5file.close()

def getLowFreqSites(args):
    seq_error = 0.0001
    #allcFiles = glob.glob(args['path_to_allc'] + "/allc_" + args['sample_id'] + "*.tsv")
    #umethfile = filter(lambda x:re.search(r"%s" % args['unMeth'], x), allcFiles)
    umethfile = args['path_to_allc'] + "/allc_" + args['sample_id'] + "_" + args['unMeth'] + ".tsv"
    if os.path.isfile(umethfile):
        conv_rate = getconv_rate_allc(umethfile)
    else:
        die("Provide a unMethylatedControl to get conversion efficiency")
    error_rate = 1 - conv_rate + seq_error
    chrs = getChrs_allc(umethfile)
    lowfreq_pval = np.zeros(0,dtype="int8")
    allcBed = []
    chrpositions = []
    for c in chrs:
        allc = args['path_to_allc'] + "/allc_" + args['sample_id'] + "_" + c + ".tsv"
        if not os.path.isfile(allc):
            die("file %s not found!" % allc)
        chrpositions.append(len(lowfreq_pval))
        log.info("reading file %s" % allc)
        bsbed = pd.read_table(allc)
        log.info("analysing %s!" % c)
        meth_inds = np.where(bsbed['methylated'] == 1)[0]
        mcpval = callMPs(np.array(bsbed['mc_count'])[meth_inds], np.array(bsbed['total'])[meth_inds], 1 - error_rate, alternative="less")
        cpval = np.zeros(bsbed.shape[0], dtype="int8")
        cpval[meth_inds[np.where(mcpval < args['pvalue_thres'])[0]]] = 1
        lowfreq_pval = np.append(lowfreq_pval, cpval)
        try:
            allcBed = pd.concat([allcBed, bsbed])
        except TypeError:
            allcBed = bsbed
        log.info("done!")
    allcBed['lowfreq'] = pd.Series(lowfreq_pval, index=allcBed.index)
    log.info("writing the data into a h5 file")
    outH5File(allcBed,chrpositions, args['outFile'])
    log.info("finished!")

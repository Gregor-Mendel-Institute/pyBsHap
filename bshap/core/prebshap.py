## Main module for pyBsHap

import pysam
from pyfaidx import Fasta
import logging
import numpy as np
import pandas as pd
import re
import h5py
import string
import os.path, sys
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
import json
from . import meth5py

log = logging.getLogger(__name__)

def die(msg):
  sys.stderr.write('Error: ' + msg + '\n')
  sys.exit(1)

def getChrs(inBam):
    chrs = np.array([x['SN'] for x in inBam.header['SQ']])
    chrslen = np.array([x['LN'] for x in inBam.header['SQ']])
    return( (chrs, chrslen) )

def get_reverse_complement(seq):
    old_chars = "ACGT"
    replace_chars = "TGCA"
    tab = string.maketrans(old_chars,replace_chars)
    return(seq.translate(tab)[::-1])

def decodeFlag(flag):
    bintag = str(int(bin(flag)[2:]))[::-1]
    notag = (12 - len(bintag)) * '0'
    return(bintag + notag)

def filterRead(binread):
    rflag = decodeFlag(binread.flag)
    filter_flag = False
    if len(binread.get_aligned_pairs()) == 0 or rflag[10] != '0' or rflag[8] != '0':
        ## removing the duplicates
        ## Filter the read if alignment is not good!
        filter_flag = True
    return(filter_flag)

def findIntersectbed(bin_bed, map_bed):
    # bin_bed = [290605, 290688]
    # map_bed = [290626,290711]
    bin_bed_arr = np.array(list(range(bin_bed[0], bin_bed[1])))
    map_bed_arr = np.array(list(range(map_bed[0], map_bed[1])))
    intersect_bed_arr = np.sort(np.intersect1d(bin_bed_arr, map_bed_arr))
    intersect_bed = [intersect_bed_arr[0], intersect_bed_arr[-1]+1]
    return(intersect_bed)

def parseContext(tair10, cid, pos, strand):
    try:
        if strand == '0':
            dnastring = tair10[cid][pos:pos+3].seq.encode('ascii').upper()
        elif strand == '1': ## make sure you can have to identify strand here
            dnastring = tair10[cid][pos-2:pos+1].seq.encode('ascii').upper()  ## Changed the position, different from forward
            dnastring = get_reverse_complement(dnastring)
    except:
        return("CNN", "CN")
    if dnastring[1].upper() == 'G':
        dna_context = ["CG",0]
    elif dnastring[2].upper() == 'G':
        dna_context = ["CHG",1]
    elif dnastring:
        dna_context = ["CHH",2]
    return(dnastring, dna_context)


def get_read_tags(binread, directional_bs_tag = 'XR', reference_bs_tag = 'XG'):
    """
    Tags correspond to which strand read is mapped to in the reference.
     bismark manual suggests (https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html#ii-bismark-alignment-step)
        1. XR corresponds to read conversion state
        2. XG corresponds to genome conversion state.
    please check the tags while using any other aligner than bismark
    """
    tag_options = ['CT', 'GA']
    assert directional_bs_tag in dict( binread.get_tags() ), "input bam file doesnt contain read tag %s" % directional_bs_tag
    assert reference_bs_tag in dict( binread.get_tags() ), "input bam file doesnt contain reference tag %s" % reference_bs_tag
    ref_tag = dict(binread.get_tags())[reference_bs_tag]
    read_tag = dict(binread.get_tags())[directional_bs_tag]
    return([ref_tag + read_tag, ref_tag, read_tag])

def getStrand_Bismark(strand):
    if strand == 'CT': ## Forward read, count the number of C's
        re_string = ['C','T','.']
    elif strand == 'GA': ## reverse
        re_string = ['G', 'A',',']
    return(re_string)

def getMatchedSeq(binread, tair10, gap_fill = '-'):
    refseq = ''
    rseq = ''
    for (query_ix, ref_ix) in binread.get_aligned_pairs():
        if query_ix is not None and ref_ix is not None:
            rseq = rseq + binread.seq[query_ix]
            refseq = refseq + tair10[binread.reference_id][ref_ix].seq
        elif query_ix is None and ref_ix is not None:
            rseq = rseq + gap_fill
            refseq = refseq + tair10[binread.reference_id][ref_ix].seq
        elif ref_ix is None and query_ix is not None:
            rseq = rseq + binread.seq[query_ix]
            refseq = refseq + gap_fill
    return(refseq, rseq)

def getHighlightedSeqs(refseq, rseq, strand):
    re_string = getStrand_Bismark(strand)
    dot_rseq = ''
    c_haps = ''
    for i,c in enumerate(refseq):
        if c.upper() != re_string[0]:  ### do not care about sites when there is no C in reference 
            dot_rseq = dot_rseq + re_string[2] 
        else:
            dot_rseq = dot_rseq + rseq[i]
            ### get only C's and concat them together
            if rseq[i].upper() == re_string[0]: ## the read also has a C
                c_haps = c_haps + '1'
            elif rseq[i].upper() == re_string[1]:
                c_haps = c_haps + '0'
    return((dot_rseq, c_haps))

def getSeqRecord(binread, tair10, ref_strand):
    ## bin_bed = ['Chr1', start, end, binLen, pointPos]
    char_add = 'N'
    re_string = getStrand_Bismark(ref_strand)
    refseq, rseq = getMatchedSeq(binread, tair10)
    dot_rseq, c_haps = getHighlightedSeqs(refseq, rseq, ref_strand)
    return((SeqRecord(Seq(dot_rseq, generic_dna), id = binread.query_name.split(' ')[0]), c_haps))

#### ==========================================
##      MHL calculation from here
##  adapted from Guo et al, 2017
#### ==========================================


def get_hap_comb_single(ecs_hap_bins):
    ncols = len(ecs_hap_bins)
    hap_comb_single = np.zeros(0, dtype=ecs_hap_bins.dtype)
    for i in range(ncols):
        for j in range(i, ncols):
            hap_comb_single = np.append(hap_comb_single, ecs_hap_bins[i:j+1])
    return(hap_comb_single)

def get_hap_comb(cs_hap_bins):
    hap_comb = np.zeros(0, dtype=cs_hap_bins.dtype)
    for i in range(len(cs_hap_bins)):
        hap_comb = np.append(hap_comb, get_hap_comb_single(cs_hap_bins[i]))
    return(hap_comb)

def calculate_mhl(cs_hap_bins):
    ## this is adapted from Guo et al. 2017 Figure 2.
    if len(cs_hap_bins) <= 1:
        return((np.nan, 0, 0))
    nplen = np.vectorize(len)
    hap_comb = get_hap_comb(cs_hap_bins)
    mhl_value = np.zeros(0, dtype=float)
    ncols = np.max(nplen(hap_comb))
    for i in range(ncols):
        input_tmp = np.where(nplen(hap_comb) == i + 1)[0]
        nmeth = len( np.where( hap_comb[input_tmp] == '1' * (i + 1) )[0]  )
        if len(input_tmp) >= len(cs_hap_bins)/1.3:
            mhl_value = np.append(mhl_value, (i + 1) * nmeth / float(len(input_tmp)))
    mhl_value_final = np.sum(mhl_value) / sum(range(ncols + 1))
    return((mhl_value_final, len(cs_hap_bins), ncols))

def get_mhl_entire_bed(inBam, tair10, required_bed, outstat = ''):
    window_size = required_bed[3]
    read_length_thres = window_size/2
    # bin_start = required_bed[1] - (required_bed[1] % window_size)
    bin_start = required_bed[1]# - (required_bed[1] % window_size)
    estimated_bins = list(range(bin_start, required_bed[2], window_size))
    progress_bins = 0
    #### iter meths in the required bed.
    # meths_bins = meths.iter_bed_windows([required_bed[0], bin_start, required_bed[2]], window_size)
    frac_mhl = pd.DataFrame( columns=["chr", "start", "end", "mhl"] )
    for bins in estimated_bins:        ## sliding windows with window_size
        progress_bins += 1
        cs_hap_bins_directional = {}
        for etag in ['CTCT', 'GACT', 'CTGA', 'GAGA']:  ### considering all the 4 types of tags
            cs_hap_bins_directional[etag] = np.zeros(0, dtype="str")
        cs_hap_bins = np.zeros(0, dtype="str")
        bin_bed = [required_bed[0], bins, bins + window_size - 1]
        n_reads = 0
        for binread in inBam.fetch(str(bin_bed[0]), bin_bed[1], bin_bed[2]):
            r_strand = get_read_tags(binread)
            if filterRead(binread): ## Filtering the dupplicated reads.
                continue
            rseq_record, cs_bin = getSeqRecord(binread, tair10, r_strand[1])
            if len(cs_bin) > 0:
                n_reads += 1
                cs_hap_bins_directional[r_strand[0]] = np.append(cs_hap_bins_directional[r_strand[0]], cs_bin)
                cs_hap_bins = np.append(cs_hap_bins, cs_bin)
        mhl_value, no_cs_hap_bins, ncols = calculate_mhl(cs_hap_bins)
        # wma_win = meths.MethylationSummaryStats(temp_meth_bins[1], 1)
        frac_mhl.loc[progress_bins, "chr"] = bin_bed[0]
        frac_mhl.loc[progress_bins, "start"] = bin_bed[1] + 1
        frac_mhl.loc[progress_bins, "end"] = bin_bed[2] + 1
        frac_mhl.loc[progress_bins, "mhl"] = mhl_value
        # frac_mhl.loc[progress_bins, "wma"] = wma_win
        # frac_mhl.loc[progress_bins, "num_cs"] = len( temp_meth_bins[1] )
        for etag in cs_hap_bins_directional.keys():
            e_mhl_value, e_cs_bins, e_ncols = calculate_mhl(cs_hap_bins_directional[etag])
            frac_mhl.loc[progress_bins, "mhl_" + etag ] = e_mhl_value
            frac_mhl.loc[progress_bins, "nreads_" + etag ] = cs_hap_bins_directional[etag].shape
        if outstat != '':  ### Write to a file, if given
            outstat.write("%s,%s,%s,%s,%s,%s\n" % (bin_bed[0], str(bin_bed[1]), str(bin_bed[2]), frac_mhl, no_cs_hap_bins))
        if progress_bins % 1000 == 0:
            log.info("ProgressMeter - %s windows in analysed, %s total" % (progress_bins, len(estimated_bins)))
    return(frac_mhl)

def potato_mhl_calc(args):
    log.info("loading the input files!")
    inBam = pysam.AlignmentFile(args['inFile'], "rb")
    if args['inhdf5'] != '':
        meths = meth5py.load_hdf5_methylation_file(args['inhdf5'])
    else:
        meths = ''
    (chrs, chrslen) = getChrs(inBam)
    tair10 = Fasta(args['fastaFile'])
    log.info("done!")
    if args['outFile'] != "STDOUT":
        outstat = open(args['outFile'], 'w')
        outstat.write("chr,start,end,mhl_stat,wma_win,nreads\n")
    else:
        outstat = ''
        print("chr,start,end,mhl_stat,wma_win,nreads")
    if args['reqRegion'] != '0,0,0':
        required_region = args['reqRegion'].split(',')
        required_bed = [required_region[0], int(required_region[1]), int(required_region[2]), args['window_size']]
        log.info("analysing region %s:%s-%s !" % (required_bed[0], required_bed[1], required_bed[2]))
        get_mhl_entire_bed(inBam, tair10, required_bed, outstat)
        log.info("finished!")
        return(0)
    for cid, clen in zip(chrs, chrslen):     ## chromosome wise
        log.info("analysing chromosome: %s" % cid)
        required_bed = [cid, 1, clen, window_size]   ### 0 is to not print reads
        get_mhl_entire_bed(inBam, tair10, required_bed, outstat)
    log.info("finished!")
    return(0)

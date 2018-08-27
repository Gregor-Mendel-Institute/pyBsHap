# Main module for bisulfite analysis
# Summary statistics
import logging
import h5py as h5
import numpy as np
import pandas as pd
import os
import os.path as op
import glob
import sys
from . import run_bedtools
import csv
import itertools

log = logging.getLogger(__name__)

from . import genome as g
genome = g.ArabidopsisGenome()

def die(msg):
  sys.stderr.write('Error: ' + msg + '\n')
  sys.exit(1)

def get_chrs(allc_id, allc_path):
    #files_present = glob.glob(allc_path + "/allc_" + allc_id + "*tsv")
    if op.isfile(allc_path + "/allc_" + allc_id + "_Chr1.tsv"):
        return(['Chr1','Chr2','Chr3','Chr4','Chr5'])
    elif op.isfile(allc_path + "/allc_" + allc_id + "_1.tsv"):
        return(['1','2','3','4','5'])
    else:
        die("Either the allc id, %s or allc file path, %s is wrong" % (allc_id, allc_path))

def read_allc_pandas_table(allc_file):
    sniffer = csv.Sniffer()
    if op.splitext(allc_file)[1] == '.gz':  ### For the new allc files
        bsbed = pd.read_table(allc_file, header = None, compression='gzip')
    else:
        allc = open(allc_file, 'rb')
        ifheader = sniffer.has_header(allc.read(4096))
        if ifheader:
            bsbed = pd.read_table(allc_file, header = 0)
            return(bsbed)
        else:
            bsbed = pd.read_table(allc_file, header = None)
    if len(bsbed.columns.values) == 7:
        bsbed.columns = np.array(['chr','pos','strand','mc_class','mc_count','total','methylated'])
    else:
        raise(NotImplementedError)
    return(bsbed)

def read_allc_files_chrwise(allc_id, allc_path):
    allcBed = []
    chrpositions = []
    tlen = 0
    sample_chrs = get_chrs(allc_id, allc_path)
    for c in sample_chrs:
        allc_file = allc_path + "/allc_" + allc_id + "_" + c + ".tsv"
        log.info("progress: reading %s!" % allc_file)
        chrpositions.append(tlen)
        bsbed = read_allc_pandas_table(allc_file)
        tlen = tlen + bsbed.shape[0]
        try:
            allcBed = pd.concat([allcBed, bsbed])
        except TypeError:
            allcBed = bsbed
    return((allcBed, chrpositions))

def generage_h5file_from_allc(allc_id, allc_path, outFile, pval_thres=0.01):
    if allc_path == 'new':
        log.info("reading the allc file")
        allcBed = read_allc_pandas_table(allc_id)
        log.info("done!")
        if len(np.unique(allcBed.iloc[0:100000,6])) == 1: ### Check if the top one million is all either 1's or 0's
            from . import bsseq
            pval = bsseq.callMPs_allcbed(allcBed)
            allcBed.iloc[np.where(pval > pval_thres)[0],6] = 0   ### changing the methylated column to 0 for non methylated sites.
        allcBed_new = pd.DataFrame(columns = allcBed.columns)
        chrpositions = np.zeros(1, dtype="int")
        log.info("sorting the bed file!")
        for ec in genome.chrs:
            allcBed_echr = allcBed.iloc[np.where(allcBed.iloc[:,0] == ec)[0], :]
            allcBed_new = allcBed_new.append(allcBed_echr.iloc[np.argsort(np.array(allcBed_echr.iloc[:,1], dtype=int)), :] , ignore_index=True)
            chrpositions = np.append(chrpositions, chrpositions[-1] + allcBed_echr.shape[0])
        log.info("done!")
        allcBed = allcBed_new
    else:
        (allcBed, chrpositions) = read_allc_files_chrwise(allc_id, allc_path)
    log.info("writing a hdf5 file")
    generate_H5File(allcBed,chrpositions, outFile)

chunk_size = 100000
def generate_H5File(allcBed, chrpositions, outFile):
    h5file = h5.File(outFile, 'w')
    num_lines = len(chrpositions)
    h5file.create_dataset('chunk_size', data=chunk_size, shape=(1,),dtype='i8')
    h5file.create_dataset('chrpositions', data=chrpositions, shape=(num_lines,),dtype='i4')
    allc_columns = ['chr', 'pos', 'strand', 'mc_class', 'mc_count', 'total', 'methylated']
    allc_dtypes = ['S8', 'int', 'S4', 'S8', 'int', 'int', 'int8']
    ## Going through all the columns
    for cols, coltype in zip(allc_columns, allc_dtypes):
        h5file.create_dataset(cols, compression="gzip", data=np.array(allcBed[cols],dtype=coltype), shape=(allcBed.shape[0],), chunks = ((chunk_size,)))
    if allcBed.shape[1] > 7:
        for i in range(7,allcBed.shape[1]):
            try:
                h5file.create_dataset(allcBed.columns[i], compression="gzip", data=np.array(allcBed[allcBed.columns[i]]), shape=(allcBed.shape[0],),chunks = ((chunk_size,)))
            except TypeError:
                h5file.create_dataset(allcBed.columns[i], compression="gzip", data=np.array(allcBed[allcBed.columns[i]],dtype="S8"), shape=(allcBed.shape[0],),chunks = ((chunk_size,)))
    h5file.close()



def load_hdf5_methylation_file(hdf5_file, bin_bed=''):
    return HDF5MethTable(hdf5_file, bin_bed)

# Try to make it as a class, learned from PyGWAS
class HDF5MethTable(object):

    def __init__(self,hdf5_file, bin_bed=''):
        self.h5file = h5.File(hdf5_file,'r')
        self.filter_pos_ix = self.get_filter_inds(bin_bed)
        self.chrpositions = np.array(self.h5file['chrpositions'])
        self.positions = self.h5file['pos']
        self.chr = self.h5file['chr']
        self.methylated = self.h5file['methylated']
        self.strand = self.h5file['strand']
        self.mc_class = self.h5file['mc_class']
        self.mc_count = self.h5file['mc_count']
        self.mc_total = self.h5file['total']


    def close_h5file(self):
        if self.h5file is not None:
            self.h5file.close()

    def get_filter_inds(self, bin_bed):
        # bin_bed = ['Chr1', 0, 100]
        if bin_bed == '':
            return(None)
        if len(bin_bed) != 3:
            die("provide the genomic position as array of length 3. ex., ['Chr1', 0, 100]")
        req_chr_ind, chr_inds = self.get_chrinds(bin_bed[0])
        if req_chr_ind is None:
            return(np.zeros(0, dtype=int))
        req_inds = chr_inds[0] + np.searchsorted(self.positions[chr_inds[0]:chr_inds[1]], [bin_bed[1],bin_bed[2]], side='right')
        return(np.arange(req_inds[0],req_inds[1]))

    def get_chrs_list(self, filter_pos_ix):
        if filter_pos_ix is None:
            filter_pos_ix = np.arange(len(self.positions))
        if len(filter_pos_ix) == 0:
            return(np.zeros(0, dtype="S8"))
        chr_list = np.zeros(len(filter_pos_ix), dtype="S8")
        filter_pos_ix = np.array(filter_pos_ix, dtype="int")
        for echr in genome.chrs:
            chrinds = self.get_chrinds(echr)
            chr_list[np.where((filter_pos_ix >= chrinds[1][0]) & (filter_pos_ix < chrinds[1][1]))[0]] = echr
        return(chr_list)

    def get_positions(self, filter_pos_ix=None):
        if filter_pos_ix is None:
            return(np.array(self.h5file['pos']))
        elif type(filter_pos_ix) is np.ndarray:
            ### provide a sorted filter_pos_ix
            rel_pos_ix = filter_pos_ix - filter_pos_ix[0] ### Make indices relative to chromosome
            return(self.h5file['pos'][filter_pos_ix[0]:filter_pos_ix[-1]+1][rel_pos_ix])
        else:
            return(self.h5file['pos'][filter_pos_ix])

    def get_bed_df(self, filter_pos_ix, full_bed=False, read_threshold=0):
        req_pos = self.get_positions(filter_pos_ix)
        if full_bed:
            conname = self.get_mc_class(filter_pos_ix)
            strand = self.get_strand(filter_pos_ix)
            permeth = self.get_permeths(filter_pos_ix, read_threshold=read_threshold)
            return(pd.DataFrame(np.column_stack((self.get_chrs_list(filter_pos_ix),req_pos, req_pos + 1, conname, permeth, strand)), columns=['chr', 'start', 'end', 'mc_class', 'permeth', 'strand']))
        return(pd.DataFrame(np.column_stack((self.get_chrs_list(filter_pos_ix),req_pos, req_pos + 1)), columns=['chr', 'start', 'end']))

    def get_methylated(self, filter_pos_ix=None):
        if filter_pos_ix is None:
            return(np.array(self.h5file['methylated']))
        elif type(filter_pos_ix) is np.ndarray:
            rel_pos_ix = filter_pos_ix - filter_pos_ix[0]
            return(self.h5file['methylated'][filter_pos_ix[0]:filter_pos_ix[-1]+1][rel_pos_ix])
        else:
            return(self.h5file['methylated'][filter_pos_ix])

    def get_strand(self, filter_pos_ix=None):
        if filter_pos_ix is None:
            return(np.array(self.h5file['strand']))
        elif type(filter_pos_ix) is np.ndarray:
            rel_pos_ix = filter_pos_ix - filter_pos_ix[0]
            return(self.h5file['strand'][filter_pos_ix[0]:filter_pos_ix[-1]+1][rel_pos_ix])
        else:
            return(self.h5file['strand'][filter_pos_ix])

    def get_mc_class(self, filter_pos_ix=None):
        if filter_pos_ix is None:
            return(np.array(self.h5file['mc_class']))
        elif type(filter_pos_ix) is np.ndarray:
            rel_pos_ix = filter_pos_ix - filter_pos_ix[0]
            return(self.h5file['mc_class'][filter_pos_ix[0]:filter_pos_ix[-1]+1][rel_pos_ix])
        else:
            return(self.h5file['mc_class'][filter_pos_ix])

    def get_req_mc_class_ix(self, req_context, filter_pos_ix):
        if req_context is None:
            if filter_pos_ix is not None:
                return(np.arange(len(filter_pos_ix)))
            else:
                return(None)
        import re
        cor = re.compile(req_context)
        np_vmatch = np.vectorize(lambda x:bool(cor.match(x)))
        return(np.where(np_vmatch(self.get_mc_class(filter_pos_ix)))[0])

    def get_mc_count(self, filter_pos_ix=None):
        if filter_pos_ix is None:
            return(np.array(self.h5file['mc_count']))
        elif type(filter_pos_ix) is np.ndarray:
            rel_pos_ix = filter_pos_ix - filter_pos_ix[0]
            return(self.h5file['mc_count'][filter_pos_ix[0]:filter_pos_ix[-1]+1][rel_pos_ix])
        else:
            return(self.h5file['mc_count'][filter_pos_ix])

    def get_lowfreq(self, filter_pos_ix=None):
        if filter_pos_ix is None:
            return(np.array(self.h5file['lowfreq']))
        elif type(filter_pos_ix) is np.ndarray:
            rel_pos_ix = filter_pos_ix - filter_pos_ix[0]
            return(self.h5file['lowfreq'][filter_pos_ix[0]:filter_pos_ix[-1]+1][rel_pos_ix])
        else:
            return(self.h5file['lowfreq'][filter_pos_ix])

    def get_mc_total(self, filter_pos_ix=None):
        if filter_pos_ix is None:
            return(np.array(self.h5file['total']))
        elif type(filter_pos_ix) is np.ndarray:
            rel_pos_ix = filter_pos_ix - filter_pos_ix[0]
            return(self.h5file['total'][filter_pos_ix[0]:filter_pos_ix[-1]+1][rel_pos_ix])
        else:
            return(self.h5file['total'][filter_pos_ix])

    def get_permeths(self, filter_pos_ix=None, read_threshold=0):
	# If read_threshold is given
	#   if methylated then we have a value
	#   else the value is 0
	# Else the notation is -1
        methylated_cs = np.array(self.get_methylated(filter_pos_ix), dtype="float")
        mc_count = self.get_mc_count(filter_pos_ix)
        mc_total = self.get_mc_total(filter_pos_ix)
        calculated_permeths = np.divide(np.multiply(methylated_cs, mc_count), mc_total)
        if read_threshold == 0:
            return(calculated_permeths)
        else:
            if filter_pos_ix is None:
                filter_pos_ix = np.arange(len(mc_total))
            permeths = np.repeat(-1.0, len(filter_pos_ix))
            new_filter_pos_ix = np.where(mc_total > read_threshold)[0]
            permeths[new_filter_pos_ix] = calculated_permeths[new_filter_pos_ix]
            return(permeths)

    def get_chrinds(self, chrid):
        chrpositions = np.append(np.array(self.h5file['chrpositions']), len(self.h5file['pos']))
        req_chr_ind = genome.get_chr_ind(chrid)
        if req_chr_ind is None:
            return( (None, None) )
        chr_inds = [chrpositions[req_chr_ind], chrpositions[req_chr_ind + 1]]
        return((req_chr_ind, chr_inds))

    def iter_chr_windows(self, chrid, window_size):
        req_chr_ind, chr_inds = self.get_chrinds(chrid)
        return(self.iter_bed_windows([chrid, 1, genome.golden_chrlen[req_chr_ind]], window_size))

    def iter_bed_windows(self, required_bed, window_size):
        ## required_bed = ["Chr1", 1, 100]
        filter_pos_ix = self.get_filter_inds(required_bed)
        if len(filter_pos_ix) > 0:
            filter_pos = self.get_positions(filter_pos_ix)
            ind = 0
            for t in range(required_bed[1], required_bed[2], window_size):
                result = []
                bin_bed = [int(t), int(t) + window_size - 1]
                for epos in filter_pos[ind:]:
                    if epos >= bin_bed[0]:
                        if epos <= bin_bed[1]:
                            result.append(filter_pos_ix[ind])
                        elif epos > bin_bed[1] :
                            yield((bin_bed, result))
                            break
                        ind = ind + 1

def generate_meths_in_windows(meths, out_file, window_size, category=1, req_context=None):
    ## Methylation category here by default is weighted average
    # look the function below (MethylationSummaryStats).
    if op.isfile(out_file):
        die("ouput bedgraph file (%s) is already present" % out_file)
    outmeths_avg = open(out_file, 'w')
    for echr, echrlen in zip(genome.chrs, genome.golden_chrlen):
        self_windows = meths.iter_chr_windows(echr, window_size)
        count = 0
        log.info("analyzing %s" % echr)
        for ewin in self_windows:
            req_meth_avg = MethylationSummaryStats(meths, ewin[1], category, req_context)
            outmeths_avg.write("%s\t%s\t%s\t%s\n" % (echr, ewin[0][0], ewin[0][1], req_meth_avg))
            count = count + 1
            if count % 1000 == 0:
                log.info("progress: analysed %s windows" % count)
    outmeths_avg.close()
    log.info("done!")

def expand_nucleotide_code(mc_type):
    iub_dict = {"N":["A","C","G","T"],
                "H":["A","C","T"],
                "D":["A","G","T"],
                "B":["C","G","T"],
                "A":["A","C","G"],
                "R":["A","G"],
                "Y":["C","T"],
                "K":["G","T"],
                "M":["A","C"],
                "S":["G","C"],
                "W":["A","T"],
                "C":["C"],
                "G":["G"],
                "T":["T"],
                "A":["A"]}

    mc_class = list(mc_type) # copy
    if "C" in mc_type:
        mc_class.extend(["CGN", "CHG", "CHH","CNN"])
    elif "CG" in mc_type:
        mc_class.extend(["CGN"])

    for motif in mc_type:
        mc_class.extend(["".join(i) for i in
                         itertools.product(*[iub_dict[nuc] for nuc in motif])])
    return(set(mc_class))

def MethylationSummaryStats(meths, filter_pos_ix, category, req_context = None):
    # meths is already filtered for bin_bed positions
    if len(filter_pos_ix) == 0:
        return(np.nan)
    # Filtering for context
    if req_context is not None:
        filter_pos_ix = np.array(filter_pos_ix, dtype=int)[meths.get_req_mc_class_ix(req_context, filter_pos_ix)]
        if len(filter_pos_ix) == 0:
            return(np.nan)
    mc_total = meths.get_mc_total(filter_pos_ix)
    mc_count = meths.get_mc_count(filter_pos_ix)
    if np.sum(mc_total) == 0:
        return(np.nan)
    if category == 1:   # weighted mean
        methylated_cs = np.where(meths.get_methylated(filter_pos_ix) == 1)[0]
        return(np.divide(float(np.sum(mc_count[methylated_cs])), np.sum(mc_total)))
    elif category == 2: # fraction of methylated positions
        methylated = meths.get_methylated(filter_pos_ix)
        meths_len = np.sum(methylated)
        return(float(meths_len)/len(filter_pos_ix))
    elif category == 3: ## absolute means
        return(np.mean(np.divide(np.array(mc_count, dtype="float"), mc_total)))
    elif category == 4:  ### weighted mean but without making correction for converion.
        return(np.divide(float(np.sum(mc_count)), np.sum(mc_total)))
    else:
        raise(NotImplementedError)

def methylation_average_required_bed(meths, required_bed, outmeths_avg, category):
    ### required_bed = ["Chr1", 1, 1000]
    filter_pos_ix = meths.get_filter_inds(required_bed)
    req_meth_avg = MethylationSummaryStats(meths, filter_pos_ix, category)
    if outmeths_avg == '':
        return(req_meth_avg)
    elif outmeths_avg is not None:
        outmeths_avg.write("%s\t%s\t%s\t%s\n" % (required_bed[0], required_bed[1], required_bed[2], req_meth_avg))
    else:
        print(("%s\t%s\t%s\t%s\n" % (required_bed[0], required_bed[1] + 1, required_bed[2], req_meth_avg)))

def potatoskin_methylation_averages(args):
    # bin_bed = Chr1,1,100
    if args['allc_path'] == "bed":
        args['outFile'] =  args['outFile'] + '.bedGraph'
        run_bedtools.get_genomewide_methylation_average(args)
        log.info("finished!")
        return(0)
    elif args['allc_path'] ==  "hdf5":
        log.info("reading hdf5 file!")
        meths = load_hdf5_methylation_file(args['inFile'])
        log.info("done!")
    else:
        if args['allc_path'] == 'new':
            outhdf5 = op.splitext(op.splitext(args['inFile'])[0])[0] + ".hdf5"
        else:
            outhdf5 = args['allc_path'] + "allc_" + args['inFile'] + ".hdf5"
        if op.isfile(outhdf5):
            log.info("loading the hdf5 file %s!" % outhdf5)
        else:
            log.info("generating hdf5 file from allc files, %s!" % outhdf5)
            generage_h5file_from_allc(args['inFile'], args['allc_path'], outhdf5)
        meths = load_hdf5_methylation_file(outhdf5)
        log.info("done!")
    if args['required_region'] == '0,0,0':
        generate_meths_in_windows(meths, args['outFile'], args['window_size'], category=args['category'], req_context=args['context'])
        return(0)
    if args['outFile'] is not None:
        outmeths_avg = open(args['outFile'], 'w')
    else:
        outmeths_avg = None
    required_region = args['required_region'].split(',')
    required_bed = [required_region[0], int(required_region[1]), int(required_region[2])]
    log.info("analysing region %s:%s-%s !" % (required_bed[0], required_bed[1], required_bed[2]))
    methylation_average_required_bed(meths, required_bed, outmeths_avg, args['category'])
    if outmeths_avg is not None:
        outmeths_avg.close()
    log.info('finished!')
    return(0)

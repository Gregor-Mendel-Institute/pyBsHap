# Main module for bisulfite analysis
# Summary statistics
import logging
import h5py as h5
import numpy as np
import pandas as pd
import os.path
import glob
import sys
from . import run_bedtools
import csv
import itertools

log = logging.getLogger(__name__)
chrs = ['Chr1','Chr2','Chr3','Chr4','Chr5']
chrslen = [34964571, 22037565, 25499034, 20862711, 31270811]

def die(msg):
  sys.stderr.write('Error: ' + msg + '\n')
  sys.exit(1)

def get_chrs(allc_id, allc_path):
    #files_present = glob.glob(allc_path + "/allc_" + allc_id + "*tsv")
    if os.path.isfile(allc_path + "/allc_" + allc_id + "_Chr1.tsv"):
        return(['Chr1','Chr2','Chr3','Chr4','Chr5'])
    elif os.path.isfile(allc_path + "/allc_" + allc_id + "_1.tsv"):
        return(['1','2','3','4','5'])
    else:
        die("Either the allc id, %s or allc file path, %s is wrong" % (allc_id, allc_path))

def read_allc_pandas_table(allc_file):
    sniffer = csv.Sniffer()
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
        raise NotImplementedError
    return(bsbed)

def generage_h5file_from_allc(allc_id, allc_path, outFile):
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
    log.info("writing a hdf5 file")
    generate_H5File(allcBed,chrpositions, outFile)

chunk_size = 100000
def generate_H5File(allcBed, chrpositions, outFile):
    h5file = h5.File(outFile, 'w')
    num_lines = len(chrpositions)
    h5file.create_dataset('chrpositions', data=chrpositions, shape=(num_lines,),dtype='i4')
    h5file.create_dataset('pos', compression="gzip", data=np.array(allcBed['pos']), shape=(allcBed.shape[0],), dtype='i4')
    h5file.create_dataset('chunk_size', data=chunk_size, shape=(1,),dtype='i4')
    allc_columns = ['chr','strand','mc_class','mc_count','total','methylated']
    ## Going through all the columns
    for cols in allc_columns:
        try:
            h5file.create_dataset(cols, compression="gzip", data=np.array(allcBed[cols]), shape=(allcBed.shape[0],), chunks = ((chunk_size,)))
        except TypeError:
            h5file.create_dataset(cols, compression="gzip", data=np.array(allcBed[cols],dtype="S8"), shape=(allcBed.shape[0],), chunks = ((chunk_size,)))
    if allcBed.shape[1] > 7:
        for i in range(7,allcBed.shape[1]):
            try:
                h5file.create_dataset(allcBed.columns[i], compression="gzip", data=np.array(allcBed[allcBed.columns[i]]), shape=(allcBed.shape[0],),chunks = ((chunk_size,)))
            except TypeError:
                h5file.create_dataset(allcBed.columns[i], compression="gzip", data=np.array(allcBed[allcBed.columns[i]],dtype="S"), shape=(allcBed.shape[0],),chunks = ((chunk_size,)))
    h5file.close()

def load_hdf5_methylation_file(hdf5_file, bin_bed=''):
    return HDF5MethTable(hdf5_file, bin_bed)

def iter_inds(t_inds):
    result = []
    for t in t_inds:
        result.append(t)
        if len(result) == chunk_size:
            yield(result)
            result = []
    if result:
        yield(result)

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
        x_pos = self.h5file['pos']
        x_chrpositions = np.append(np.array(self.h5file['chrpositions']), len(x_pos))
        req_chr_ind = np.where(np.array(chrs) == bin_bed[0])[0][0]
        req_chr_pos_inds = [x_chrpositions[req_chr_ind],x_chrpositions[req_chr_ind + 1]]
        req_inds = x_chrpositions[req_chr_ind] + np.searchsorted(x_pos[req_chr_pos_inds[0]:req_chr_pos_inds[1]],[bin_bed[1],bin_bed[2]], side='right')
        return(range(req_inds[0],req_inds[1]))


    def get_chrs(self, filter_pos_ix):
        if filter_pos_ix is not None:
            inds = iter_inds(filter_pos_ix)
            for i in range(0, len(filter_pos_ix), chunk_size):
                yield(np.array(self.h5file['chr'][next(inds)]))
        else:
            yield(self.h5file['chr'])

    def get_positions(self, filter_pos_ix):
        if filter_pos_ix is not None:
            inds = iter_inds(filter_pos_ix)
            for i in range(0, len(filter_pos_ix), chunk_size):
                yield(np.array(self.h5file['pos'][next(inds)]))
        else:
            yield(self.h5file['pos'])

    def get_methylated(self, filter_pos_ix):
        if filter_pos_ix is not None:
            inds = iter_inds(filter_pos_ix)
            for i in range(0, len(filter_pos_ix), chunk_size):
                yield(np.array(self.h5file['methylated'][next(inds)]))
        else:
            yield(self.h5file['methylated'])

    def get_strand(self, filter_pos_ix):
        if filter_pos_ix is not None:
            inds = iter_inds(filter_pos_ix)
            for i in range(0, len(filter_pos_ix), chunk_size):
                yield(np.array(self.h5file['strand'][next(inds)]))
        else:
            yield(self.h5file['strand'])

    def get_mc_class(self, filter_pos_ix):
        if filter_pos_ix is not None:
            inds = iter_inds(filter_pos_ix)
            for i in range(0, len(filter_pos_ix), chunk_size):
                yield(np.array(self.h5file['mc_class'][next(inds)]))
        else:
            yield(self.h5file['mc_class'])

    def get_mc_count(self, filter_pos_ix):
        if filter_pos_ix is not None:
            inds = iter_inds(filter_pos_ix)
            for i in range(0, len(filter_pos_ix), chunk_size):
                yield(np.array(self.h5file['mc_count'][next(inds)]))
        else:
            yield(self.h5file['mc_count'])

    def get_lowfreq(self, filter_pos_ix):
        if filter_pos_ix is not None:
            inds = iter_inds(filter_pos_ix)
            for i in range(0, len(filter_pos_ix), chunk_size):
                yield(np.array(self.h5file['lowfreq'][next(inds)]))
        else:
            yield(self.h5file['lowfreq'])

    def get_mc_total(self, filter_pos_ix):
        if filter_pos_ix is not None:
            inds = iter_inds(filter_pos_ix)
            for i in range(0, len(filter_pos_ix), chunk_size):
                yield(np.array(self.h5file['total'][next(inds)]))
        else:
            yield(self.h5file['total'])


def MethylationSummaryStats(meths, filter_pos_ix, category):
    # meths is already filtered for bin_bed positions
    mc_total = meths.mc_total[filter_pos_ix]
    mc_count = meths.mc_count[filter_pos_ix]
    if len(filter_pos_ix) == 0:
        return(None)
    if np.sum(mc_total) == 0:
        return(None)
    if category == 1:   # weighted mean
        return np.divide(float(np.sum(mc_count)), np.sum(mc_total))
    elif category == 2: # fraction of methylated positions
        methylated = meths.methylated[filter_pos_ix]
        meths_len = np.sum(methylated)
        return float(meths_len)/len(filter_pos_ix)
    elif category == 3: ## absolute means
        return(np.mean(np.divide(np.array(mc_count, dtype="float"), mc_total)))
    else:
        raise NotImplementedError

def get_Methlation_required_bed(meths, required_bed, binLen, outmeths_avg, category):
    bin_start = required_bed[1] - (required_bed[1] % binLen)
    estimated_bins = range(bin_start, required_bed[2], binLen)
    log.info("writing methylation summary stats for %s windows!" % len(estimated_bins))
    bin_count = 0
    for bins in estimated_bins:
        bin_count = bin_count + 1
        bin_bed = [required_bed[0], bins, bins + binLen]
        filter_pos_ix = meths.get_filter_inds(bin_bed)
        if len(filter_pos_ix) > 0:
            req_meth_avg = MethylationSummaryStats(meths, filter_pos_ix, category)
        if outmeths_avg is not None:
            outmeths_avg.write("%s\t%s\t%s\t%s\n" % (bin_bed[0], bin_bed[1] + 1, bin_bed[2], req_meth_avg))
        else:
            print("%s\t%s\t%s\t%s\n" % (bin_bed[0], bin_bed[1] + 1, bin_bed[2], req_meth_avg))
        if bin_count % 100 == 0:
            log.info("progress: %s windows" % bin_count)
    return 0

def get_Methlation_GenomicRegion(args):
    # bin_bed = Chr1,1,100
    if args['required_region'] == '0,0,0':
        outFile =  args['outFile'] + '.bedGraph'
        if args['window_size'] is None:
            window_size = 200
        else:
            window_size = int(args['window_size'])
        run_bedtools.get_genomewide_methylation_WeightedMean(args['bedtoolsPath'], args['inFile'], outFile, window_size=window_size, overlap=args['overlap'], category = args['category'])
        log.info("finished!")
        return 0
    if args['allc_path'] is None:
        log.info("reading hdf5 file!")
        meths = load_hdf5_methylation_file(args['inFile'])
        log.info("done!")
    else:
        outhdf5 = args['allc_path'] + "allc_" + args['inFile'] + ".hdf5"
        if os.path.isfile(outhdf5):
            log.info("loading the hdf5 file %s!" % outhdf5)
        else:
            log.info("generating hdf5 file from allc files, %s!" % outhdf5)
            generage_h5file_from_allc(args['inFile'], args['allc_path'], outhdf5)
        meths = load_hdf5_methylation_file(outhdf5)
        log.info("done!")
    if args['outFile'] is not None:
        outmeths_avg = open(args['outFile'], 'w')
    else:
        outmeths_avg = None
    required_region = args['required_region'].split(',')
    required_bed = [required_region[0], int(required_region[1]), int(required_region[2])]
    log.info("analysing region %s:%s-%s !" % (required_bed[0], required_bed[1], required_bed[2]))
    if args['window_size'] is not None:
        get_Methlation_required_bed(meths, required_bed, args['window_size'], outmeths_avg, args['category'])
    else:
        binLen = int(required_region[2]) - int(required_region[1]) + 1
        get_Methlation_required_bed(meths, required_bed, binLen, outmeths_avg, args['category'])
    if outmeths_avg is not None:
        outmeths_avg.close()
    log.info('finished!')
    return(0)

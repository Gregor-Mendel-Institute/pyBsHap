# Main module for bisulfite analysis
# Summary statistics
import logging
import h5py as h5
import numpy as np
import pandas as pd
import os.path

log = logging.getLogger(__name__)
chrs = ['Chr1','Chr2','Chr3','Chr4','Chr5']
chrslen = [34964571, 22037565, 25499034, 20862711, 31270811]

def generage_h5file_from_allc(allc_id, allc_path, outFile):
    allcBed = []
    chrpositions = []
    tlen = 0
    for c in chrs:
        allc = allc_path + "/allc_" + allc_id + "_" + c + ".tsv"
        log.info("progress: reading %s!" % allc)
        bsbed = pd.read_table(allc)
        tlen = tlen + bsbed.shape[0]
        chrpositions.append(tlen)
        try:
            allcBed = pd.concat([allcBed, bsbed])
        except TypeError:
            allcBed = bsbed
    log.info("writing a h5 file")
    generate_H5File(allcBed,chrpositions, outFile)

def generate_H5File(allcBed, chrpositions, outFile):
    h5file = h5.File(outFile, 'w')
    num_lines = len(chrpositions)
    chunk_size = 10000
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
        x_chrpositions = np.array(self.h5file['chrpositions'])
        x_pos = self.h5file['pos']
        req_chr_ind = np.where(np.array(chrs) == bin_bed[0])[0][0]
        if req_chr_ind == 0:
            req_chr_pos_inds = [0,x_chrpositions[req_chr_ind]]
        else:
            req_chr_pos_inds = [x_chrpositions[req_chr_ind-1],x_chrpositions[req_chr_ind]]
        req_inds = np.searchsorted(x_pos[req_chr_pos_inds[0]:req_chr_pos_inds[1]],[bin_bed[1],bin_bed[2]], side='right')
        return(range(req_inds[0],req_inds[1]))


    def get_chrs(self, filter_pos_ix):
        if filter_pos_ix is not None:
            return(np.array(self.h5file['chr'][filter_pos_ix]))
        else:
            return(self.h5file['chr'])

    def get_positions(self, filter_pos_ix):
        if filter_pos_ix is not None:
            return(np.array(self.h5file['pos'][filter_pos_ix]))
        return(self.h5file['pos'])

    def get_methylated(self, filter_pos_ix):
        if filter_pos_ix is not None:
            return(np.array(self.h5file['methylated'][filter_pos_ix]))
        else:
            return(self.h5file['methylated'])

    def get_strand(self, filter_pos_ix):
        if filter_pos_ix is not None:
            return(np.array(self.h5file['strand'][filter_pos_ix]))
        else:
            return(self.h5file['strand'])

    def get_mc_class(self, filter_pos_ix):
        if filter_pos_ix is not None:
            return(np.array(self.h5file['mc_class'][filter_pos_ix]))
        else:
            return(self.h5file['mc_class'])

    def get_mc_count(self, filter_pos_ix):
        if filter_pos_ix is not None:
            return(np.array(self.h5file['mc_count'][filter_pos_ix]))
        else:
            return(self.h5file['mc_count'])

    def get_lowfreq(self, filter_pos_ix):
        if filter_pos_ix is not None:
            return(np.array(self.h5file['lowfreq'][filter_pos_ix]))
        else:
            return(self.h5file['lowfreq'])

    def get_mc_total(self, filter_pos_ix):
        if filter_pos_ix is not None:
            return(np.array(self.h5file['total'][filter_pos_ix]))
        else:
            return(self.h5file['total'])

def MethylationSummaryStats(meths, category):
    # meths is already filtered for bin_bed positions
    try:
        bin_len = len(meths.__dict__['mc_total'])
    except KeyError:
        meths.get_mc_total()
        meths.get_mc_count()
    if category == 1:   # weighted mean
        return np.divide(float(np.sum(meths.mc_count)), np.sum(meths.mc_total))
    elif category == 2: # get only methylated positions
        try:
            meths_len = np.sum(meths.__dict__['methylated'])
        except KeyError:
            meths.get_methylated()
            meths_len = np.sum(meths.__dict__['methylated'])
        return float(meths_len)/bin_len
    elif category == 3: ## absolute means
        return(np.mean(np.divide(meths.mc_count, meths.mc_total)))

def get_Methlation_required_bed(meths, required_bed, binLen, outmeths_avg, category = 1):
    bin_start = required_bed[1] - (required_bed[1] % binLen)
    estimated_bins = range(bin_start, required_bed[2], binLen)
    log.info("writing methylation summary stats per window!")
    for bins in estimated_bins:
        bin_bed = [required_bed[0], bins, bins + binLen]
        meths.filter_pos_ix = meths.get_filter_inds(bin_bed)
        import ipdb; ipdb.set_trace()
        req_meth_avg = MethylationSummaryStats(meths, category)
        outmeths_avg.write("%s,%s,%s,%s\n" % (bin_bed[0], bin_bed[1], bin_bed[2], req_meth_avg))
    return 0

def get_Methlation_GenomicRegion(args):
    # bin_bed = Chr1,1,100
    binLen = int(args['window_size'])
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
    outmeths_avg = open('meths.' + args['outFile'] + '.summary.txt', 'w')
    if args['required_region'] == '0,0,0':
        for cid, clen in zip(chrs, chrslen):     ## chromosome wise
            log.info("analysing chromosome: %s" % cid)
            required_bed = [cid, 0, clen, binLen]   ### 0 is to not print reads
            get_Methlation_required_bed(meths, required_bed, binLen, outmeths_avg)
            log.info("finished!")
            return 0
    else:
        required_region = args['required_region'].split(',')
        required_bed = [required_region[0], int(required_region[1]), int(required_region[2])]
        log.info("analysing region %s:%s-%s !" % (required_bed[0], required_bed[1], required_bed[2]))
        get_Methlation_required_bed(meths, required_bed, binLen, outmeths_avg)
    log.info('finished!')
    return(0)

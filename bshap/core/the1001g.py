# Main module for analysing 1010 genomes data
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
golden_chrlen = [30427671, 19698289, 23459830, 18585056, 26975502]

def die(msg):
  sys.stderr.write('Error: ' + msg + '\n')
  sys.exit(1)

def read_bg_file(bg_file):
    sniffer = csv.Sniffer()
    t_ef = open(bg_file, 'rb')
    ifheader = sniffer.has_header(t_ef.read(4096))
    t_ef.next() #### need to do this so that we get entire line in the next read
    delimiter = sniffer.sniff(t_ef.next()).delimiter
    if ifheader:
        e_bg = pd.read_table(bg_file, header = 0, sep= delimiter)
    else:
        e_bg = pd.read_table(bg_file, header = None, sep= delimiter)
    e_bg.columns = np.array(['chr','start','end','value'])
    return(e_bg)

chunk_size = 1000
def generate_h5_1001g(bg_files, outH5file):
    ## Make sure you have exactly same number of lines in all the files
    num_lines = len(bg_files)
    log.info("writing a h5 file for %s bdg files into %s" % (num_lines, outH5file))
    h5file = h5.File(outH5file, 'w')
    h5file.create_dataset('chunk_size', data=chunk_size, shape=(1,),dtype='i8')
    h5file.create_dataset('files', data=np.array(bg_files), shape=(num_lines,))
    h5file.create_dataset('accessions', shape= ((num_lines,)), dtype='S12')
    base_bg = read_bg_file(bg_files[0])
    h5file['accessions'][0] = os.path.basename(bg_files[0]).split(".")[1]
    n_rows = base_bg.shape[0]
    h5file.create_dataset('value', shape = ((n_rows, num_lines)), chunks=((chunk_size,num_lines)),dtype='float')
    for ef_ind in range(num_lines - 1):
        e_bg = read_bg_file(bg_files[ef_ind + 1])
        h5file['accessions'][ef_ind + 1] = os.path.basename(bg_files[ef_ind + 1]).split(".")[1] ## generally accession is second in the file name, eg., mhl.10001.txt
        if n_rows == 0:
            n_rows = e_bg.shape[0]
        else:
            if n_rows != e_bg.shape[0]:
                die("please provide bedgraph files with exactly same number of lines")
        h5file['value'][:,ef_ind+1] = np.array(e_bg['value'], dtype='float')
        if ef_ind % 50 == 0 and ef_ind > 0:
            log.info("written %s files into h5file" % ef_ind)
    h5file.create_dataset('chr', shape=(n_rows,), data = np.array(base_bg['chr'], dtype='str'))
    h5file.create_dataset('start', shape=(n_rows,), data = np.array(base_bg['start'], dtype='int'))
    h5file.create_dataset('end', shape=(n_rows,), data = np.array(base_bg['end'], dtype='int'))
    h5file['value'][:,0] = np.array(base_bg['value'])
    h5file.close()
    log.info("done")

### Class object for the the1001g matrices

def load_the1001g_hdf5_value_file(hdf5_file):
    return HDF51001gTable(hdf5_file)

# Try to make it as a class, learned from PyGWAS
class HDF51001gTable(object):

    def __init__(self,hdf5_file):
        self.h5file = h5.File(hdf5_file,'r')
        self.chr = self.h5file['chr']
        self.start = self.h5file['start']
        self.end = self.h5file['end']
        self.value = self.h5file['value']
        self.accessions = np.array(self.h5file['accessions'])

    def get_chr(self, filter_pos_ix=None):
        if filter_pos_ix is None:
            return(np.array(self.h5file['chr']))
        elif type(filter_pos_ix) is np.ndarray:
            rel_pos_ix = filter_pos_ix - filter_pos_ix[0]
            return(self.h5file['chr'][filter_pos_ix[0]:filter_pos_ix[-1]+1][rel_pos_ix])
        else:
            return(self.h5file['chr'][filter_pos_ix])

    def get_start(self, filter_pos_ix=None):
        if filter_pos_ix is None:
            return(np.array(self.h5file['start']))
        elif type(filter_pos_ix) is np.ndarray:
            rel_pos_ix = filter_pos_ix - filter_pos_ix[0]
            return(self.h5file['start'][filter_pos_ix[0]:filter_pos_ix[-1]+1][rel_pos_ix])
        else:
            return(self.h5file['start'][filter_pos_ix])

    def get_end(self, filter_pos_ix=None):
        if filter_pos_ix is None:
            return(np.array(self.h5file['end']))
        elif type(filter_pos_ix) is np.ndarray:
            rel_pos_ix = filter_pos_ix - filter_pos_ix[0]
            return(self.h5file['end'][filter_pos_ix[0]:filter_pos_ix[-1]+1][rel_pos_ix])
        else:
            return(self.h5file['end'][filter_pos_ix])

    def get_value(self, filter_pos_ix=None):
        if filter_pos_ix is None:
            return(np.array(self.h5file['value']))
        elif type(filter_pos_ix) is np.ndarray:
            rel_pos_ix = filter_pos_ix - filter_pos_ix[0]
            return(self.h5file['value'][filter_pos_ix[0]:filter_pos_ix[-1]+1,:][rel_pos_ix,:])
        else:
            return(self.h5file['value'][filter_pos_ix,:])

    def get_bed_df(self, filter_pos_ix):
        req_chr = self.get_chr(filter_pos_ix)
        req_start = self.get_start(filter_pos_ix)
        req_end = self.get_end(filter_pos_ix)
        return(pd.DataFrame(np.column_stack((req_chr,req_start, req_end)), columns=['chr', 'start', 'end']))

    def get_inds_overlap_region(self, region_bed):
        ## region_bed = ['Chr1',3631, 5899]
        # returns indices for the overlap of given region
        # returns even if overlap is 1 bp
        # also make sure start is always greater than end
        chr_inds = np.where(self.get_chr(None) == region_bed[0])[0]
        #chr_start = self.get_start(chr_inds)
        chr_end = self.get_end(chr_inds)
        return(np.where((chr_end >= region_bed[1]) & (chr_end < region_bed[2]))[0])

    def get_inds_overlap_region_file(self, region_file, just_names=False, araport11_file=None):
        whole_bed = self.get_bed_df(None)
        return(run_bedtools.get_filter_bed_ix(region_file, whole_bed, just_names=just_names, araport11_file=araport11_file))

    def get_inds_matching_region(self, region_bed):
        # returns values given region_bed (a pd dataframe with chr, start and end)
        ## maybe useful for wma hdf5 files
        chr_inds = np.where(self.get_chr(None) == region_bed[0])[0]
        chr_start = self.get_start(chr_inds)
        chr_end = self.get_end(chr_inds)
        return(np.where( (chr_start == region_bed[1]) & (chr_end == region_bed[2]) )[0])

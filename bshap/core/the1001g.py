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
import pybedtools as pybed

class ArabidopsisGenome(object):
    ## coordinates for ArabidopsisGenome using TAIR 10

    def __init__(self):
        self.chrs = ['Chr1','Chr2','Chr3','Chr4','Chr5']
        self.real_chrlen = [34964571, 22037565, 25499034, 20862711, 31270811]
        self.golden_chrlen = [30427671, 19698289, 23459830, 18585056, 26975502]
        self.chr_inds = np.append(0, np.cumsum(self.golden_chrlen))
        self.centro_start = [14364752, 3602775, 12674550, 2919690, 11668616]
        self.centro_end   = [15750321, 3735247, 13674767, 4011692, 12082583]
        self.cetro_mid = np.add(self.centro_start, self.centro_end)/2

    def get_bed_ids_str(self, **kwargs):
        for req_name in kwargs:
            req_bed_df = pd.read_table( kwargs[req_name], header=None )
            setattr(self, req_name, req_bed_df)
            setattr(self, req_name + "_str", np.array(req_bed_df.iloc[:,0] + ',' + req_bed_df.iloc[:,1].map(str) + ',' + req_bed_df.iloc[:,2].map(str), dtype="str") )

    def get_chr_ind(self, echr):
        echr_num = str(echr).replace("Chr", "").replace("chr", "")
        real_chrs = np.array( [ ec.replace("Chr", "").replace("chr", "") for ec in self.chrs ] )
        try:
            return(np.where(real_chrs == echr_num)[0][0])
        except IndexError, err_idx:
            return(None)


log = logging.getLogger(__name__)

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
    h5file.create_dataset('accessions', shape= ((num_lines,)), dtype='S20')
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

    def get_bed_df(self, filter_pos_ix, return_str = False):
        req_chr = self.get_chr(filter_pos_ix)
        req_start = self.get_start(filter_pos_ix)
        req_end = self.get_end(filter_pos_ix)
        if return_str:
            return(np.array(pd.Series(req_chr).map(str) + ',' + pd.Series(req_start).map(str) + ',' + pd.Series(req_end).map(str), dtype = "str" ))
        return(pd.DataFrame(np.column_stack((req_chr,req_start, req_end)), columns=['chr', 'start', 'end']))

    def get_inds_overlap_region(self, region_bed, g = None):
        ## region_bed = ['Chr1',3631, 5899]
        region_bedpy = pybed.BedTool('%s %s %s' % (region_bed[0], region_bed[1], region_bed[2]), from_string=True)
        chr_inds = np.where(self.get_chr(None) == region_bed[0])[0]
        chr_df = self.get_bed_df(chr_inds)
        if g is None:
            chr_intersect_df = pybed.BedTool.from_dataframe(chr_df).intersect(region_bedpy, wa=True).to_dataframe()
        else:
            chr_intersect_df = pybed.BedTool.from_dataframe(chr_df).intersect(region_bedpy, wa=True, sorted = True, g = g ).to_dataframe()
        chr_intersect_str = np.array(chr_intersect_df.iloc[:,0] + "," + chr_intersect_df.iloc[:,1].map(str) + "," +  chr_intersect_df.iloc[:,2].map(str), dtype="str")
        chr_str = np.array(chr_df.iloc[:,0].map(str) + ',' + chr_df.iloc[:,1].map(str) + ',' + chr_df.iloc[:,2].map(str), dtype = "str" )
        return(np.where( np.in1d( chr_str, chr_intersect_str ) )[0])

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

    def get_inds_matching_region_file(self, region_file):
        region_df = pd.read_table(region_file)
        region_cols = np.array(region_df.iloc[:,0] + "," + region_df.iloc[:,1].map(str) + "," +  region_df.iloc[:,2].map(str), dtype="str")
        bed_str = self.get_bed_df(None, return_str=True)
        return(np.where( np.in1d(bed_cols, region_cols)  )[0])

    def get_phenos_df(self, filter_pos_ix, outfile):
        values_df = np.nanmean( self.get_value(filter_pos_ix), axis = 0 )
        out = open(outfile, 'w')
        out.write("accessionid,pheno\n")
        for a_ind in range(len(self.accessions)):
            if ~np.isnan(values_df[a_ind]):
                out.write("%s,%s\n" % ( self.accessions[a_ind], values_df[a_ind] ))
        out.close()

    def matchAcc_ix(self, accs, return_np=False):
        acc_ix = [ np.where(self.accessions == i)[0][0] for i in accs]
        if return_np:
            return(np.array(acc_ix))
        return(acc_ix)

    @staticmethod
    def getInds_araGenome(df_str):
        tair10 = ArabidopsisGenome()
        if type(df_str) is not pd.core.series.Series and type(df_str) is not pd.core.frame.DataFrame:
            die("please input pandas dataframe or series object")
        elif type(df_str) is pd.core.series.Series:
            df_str_np = np.array(df_str, dtype="string")
            df_str_unique = np.unique(df_str_np, return_inverse=True)
            df_str_inds = np.array(pd.Series(df_str_unique[0]).str.split(",").apply(getInd_bin_bed, args= (tair10,)  ))
            return( df_str_inds[df_str_unique[1]] )
        elif type(df_str) is pd.core.frame.DataFrame:
            ## here first column is chr and second is position
            if df_str.shape[1] == 3:
                df_str = pd.DataFrame(df_str.iloc[:,0]).join(pd.DataFrame( ((df_str.iloc[:,1] + df_str.iloc[:,2]) / 2).apply(int) ))
            chrom = np.char.replace(np.core.defchararray.lower(np.array(df_str.iloc[:,0], dtype="string")), "chr", "")
            return(tair10.chr_inds[np.array(chrom, dtype=int) - 1] + np.array(df_str.iloc[:,1]) )

def getInd_bin_bed(bin_bed, tair10):
    ## bin_bed = ["Chr1", 1, 1000]
    bin_s = [int(bin_bed[0].replace("Chr", "")) - 1, int(bin_bed[1]), int(bin_bed[2])]
    return(tair10.chr_inds[bin_s[0]] + int(( bin_s[1] + bin_s[2] )/2) )

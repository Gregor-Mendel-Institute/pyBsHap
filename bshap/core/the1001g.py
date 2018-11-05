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

from . import genome as g
genome = g.ArabidopsisGenome()


log = logging.getLogger(__name__)

def die(msg):
  sys.stderr.write('Error: ' + msg + '\n')
  sys.exit(1)

def read_bg_file(bg_file):
    sniffer = csv.Sniffer()
    t_ef = open(bg_file, 'rb')
    ifheader = sniffer.has_header(t_ef.read(4096))
    next(t_ef) #### need to do this so that we get entire line in the next read
    delimiter = sniffer.sniff(next(t_ef)).delimiter
    if ifheader:
        e_bg = pd.read_table(bg_file, header = 0, sep= delimiter)
    else:
        e_bg = pd.read_table(bg_file, header = None, sep= delimiter)
    e_bg.columns = np.array(['chr','start','end','value'])
    return(e_bg)

def generate_h5_1001g(bg_files, outH5file, chunk_size=1000):
    ## Make sure you have exactly same number of lines in all the files
    num_lines = len(bg_files)
    log.info("writing a h5 file for %s bdg files into %s" % (num_lines, outH5file))
    h5file = h5.File(outH5file, 'w')
    h5file.create_dataset('files', data=np.array(bg_files), shape=(num_lines,))
    h5file.create_dataset('accessions', shape= ((num_lines,)), dtype='S20')
    base_bg = read_bg_file(bg_files[0])
    h5file['accessions'][0] = os.path.basename(bg_files[0]).split(".")[1]
    n_rows = base_bg.shape[0]
    if n_rows < chunk_size:
        chunk_size = n_rows
    h5file.create_dataset('chunk_size', data=chunk_size, shape=(1,),dtype='i8')
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
        self.chunk_size = self.h5file['chunk_size'][0]
        self.accessions = np.array(self.h5file['accessions'])

    def __getattr__(self, name, filter_pos_ix=None, return_np=False):
        if name not in ['chr', 'start', 'end', 'value']:
            raise AttributeError("%s is not in the keys for HDF5. Only accepted values are ['chr', 'start', 'end', 'value']" % name)
        if filter_pos_ix is None:
            if return_np:
                return(np.array(self.h5file[str(name)]))
            return(self.h5file[str(name)])
        elif type(filter_pos_ix) is np.ndarray:
            rel_pos_ix = filter_pos_ix - filter_pos_ix[0]
            ret_attr = np.array(self.h5file[str(name)][filter_pos_ix[0]:filter_pos_ix[-1]+1][rel_pos_ix])
        else:
            ret_attr = np.array(self.h5file[str(name)][filter_pos_ix])
        if name in ['chr']:
            ret_attr = ret_attr.astype('U13')
        elif name in ['start', 'end']:
            ret_attr = ret_attr.astype(int)
        else:
            ret_attr = ret_attr.astype(float)
        return(ret_attr)

    def change_accs_ids(self, sra_table):
        ## /projects/cegs/rahul/009.1001methylomes.rawdata/final_dataset_1001meths_rawdata.txt
        sra_table_1001g = pd.read_table(sra_table, header = None)
        self.accessions = np.array([ str(sra_table_1001g.iloc[:,0][np.where(sra_table_1001g.iloc[:,1] == es)[0][0]]) for es in self.accessions])

    def get_bed_df(self, filter_pos_ix, return_str = False):
        req_chr = self.__getattr__('chr', filter_pos_ix, return_np=True)
        req_start = self.__getattr__('start', filter_pos_ix, return_np=True)
        req_end = self.__getattr__('end', filter_pos_ix, return_np=True)
        if return_str:
            return(np.array(pd.Series(req_chr).map(str) + ',' + pd.Series(req_start).map(str) + ',' + pd.Series(req_end).map(str), dtype = "str" ))
        return(pd.DataFrame(np.column_stack((req_chr,req_start, req_end)), columns=['chr', 'start', 'end']))

    def get_inds_overlap_region(self, region_bed, g = None):
        ## region_bed = ['Chr1',3631, 5899]
        ## or region_bed = "Chr1,3631,5899"
        if isinstance(region_bed, basestring):
            t_split = region_bed.split(",")
            assert len(t_split) == 3
            region_bed = [ t_split[0], int(t_split[1]), int(t_split[2]) ]
        region_bedpy = pybed.BedTool('%s %s %s' % (region_bed[0], region_bed[1], region_bed[2]), from_string=True)
        chr_inds = np.where(self.__getattr__('chr', return_np=True) == region_bed[0])[0]
        chr_df = self.get_bed_df(chr_inds)
        if g is None:
            chr_intersect_df = pybed.BedTool.from_dataframe(chr_df).intersect(region_bedpy, wa=True).to_dataframe()
        else:
            chr_intersect_df = pybed.BedTool.from_dataframe(chr_df).intersect(region_bedpy, wa=True, sorted = True, g = g ).to_dataframe()
        chr_intersect_str = np.array(chr_intersect_df.iloc[:,0] + "," + chr_intersect_df.iloc[:,1].map(str) + "," +  chr_intersect_df.iloc[:,2].map(str), dtype="str")
        chr_str = np.array(chr_df.iloc[:,0].map(str) + ',' + chr_df.iloc[:,1].map(str) + ',' + chr_df.iloc[:,2].map(str), dtype = "str" )
        return(np.where( np.in1d( chr_str, chr_intersect_str ) )[0] + chr_inds[0])

    def get_inds_overlap_region_file(self, region_file, just_names=False, araport11_file=None):
        whole_bed = self.get_bed_df(None)
        return(run_bedtools.get_filter_bed_ix(region_file, whole_bed, just_names=just_names, araport11_file=araport11_file))

    def get_inds_matching_region(self, region_bed):
        # returns values given region_bed (a pd dataframe with chr, start and end)
        ## maybe useful for wma hdf5 files
        ## region_bed = ['Chr1',3631, 5899]
        ## or region_bed = "Chr1,3631,5899"
        if isinstance(region_bed, basestring):
            t_split = region_bed.split(",")
            assert len(t_split) == 3
            region_bed = [ t_split[0], int(t_split[1]), int(t_split[2]) ]
        chr_inds = np.where(np.char.decode(np.array(self.__getattr__("chr", None))) == region_bed[0])[0]
        chr_start = self.__getattr__("start", chr_inds)
        chr_end = self.__getattr__("end", chr_inds)
        return(np.where( (chr_start == region_bed[1]) & (chr_end == region_bed[2]) )[0])

    def get_inds_matching_region_file(self, region_file):
        region_df = pd.read_table(region_file)
        region_cols = np.array(region_df.iloc[:,0] + "," + region_df.iloc[:,1].map(str) + "," +  region_df.iloc[:,2].map(str), dtype="str")
        bed_str = self.get_bed_df(None, return_str=True)
        return(np.where( np.in1d(bed_cols, region_cols)  )[0])

    def get_phenos_df(self, filter_pos_ix, outfile):
        values_df = np.nanmean( self.__getattr__("value", filter_pos_ix), axis = 0 )
        if outfile is None:
            values_df_pd = pd.DataFrame( np.column_stack((self.accessions, values_df)), columns = ["accessionid", "pheno"] ).dropna()
            values_df_pd['pheno'] = pd.to_numeric(values_df_pd['pheno'])
            return(values_df_pd)
        out = open(outfile, 'w')
        out.write("accessionid,pheno\n")
        for a_ind in range(len(self.accessions)):
            if ~np.isnan(values_df[a_ind]):
                out.write("%s,%s\n" % ( self.accessions[a_ind], values_df[a_ind] ))
        out.close()

    def get_matching_accs_ix(self, accs, return_np=False):
        return(matching_accessions_ix(self.accessions, accs, return_np))


def matching_accessions_ix(target_accs, accs, return_np=False):
    acc_ix = []
    for ea in accs:
        t_ix = np.where(target_accs == ea)[0]
        if len(t_ix) == 0:
            acc_ix.append(None)
        else:
            acc_ix.append(t_ix[0])
    if return_np:
        acc_ix = np.array(acc_ix)[np.where(np.not_equal(acc_ix, None))[0]].astype("int")
    return(acc_ix)


class ContextsHDF51001gTable(object):
    ### A object for all the hdf5 file for contexts (CG, CHG and CHH)

    def __init__(self, wma_path):
        self.wma_path = wma_path
        self.load_files()

    def load_files(self):
        from glob import glob
        self.cg = HDF51001gTable(glob(self.wma_path + "/" + "*.CG.hdf5")[0])
        self.chg = HDF51001gTable(glob(self.wma_path + "/" + "*.CHG.hdf5")[0])
        self.chh = HDF51001gTable(glob(self.wma_path + "/" + "*.CHH.hdf5")[0])
        self.n_cg = HDF51001gTable(glob(self.wma_path + "/" + "*.CG.count*.hdf5")[0])
        self.n_chg = HDF51001gTable(glob(self.wma_path + "/" + "*.CHG.count*.hdf5")[0])
        self.n_chh = HDF51001gTable(glob(self.wma_path + "/" + "*.CHH.count*.hdf5")[0])
        if len( glob(self.wma_path + "/" + "*.CN.hdf5") ) > 0:
            self.cn = HDF51001gTable(glob(self.wma_path + "/" + "*.CN.hdf5")[0])

    def get_filter_inds(self, req_genes_str):
        ## req_genes_str is a list of all the strings
        all_strs = self.cg.get_bed_df(filter_pos_ix=None, return_str = True)
        return(np.where( np.in1d(all_strs, req_genes_str ) )[0])

    def get_cg_chg_chh_meths(self, filter_ind):
        ## Given list of indices, function outputs average of methylations in three contexts
        if filter_ind is None or len(filter_ind) == 0:
            return((None, None))
        t_req_gene = pd.DataFrame( np.column_stack(( np.nanmean(self.cg.__getattr__('value', filter_ind), axis = 0), np.nanmean(self.chg.__getattr__('value', filter_ind), axis = 0), np.nanmean(self.chh.__getattr__('value', filter_ind), axis = 0))), columns=["CG","CHG","CHH"] )
        t_n_req_gene = pd.DataFrame( np.column_stack(( np.nanmean(self.n_cg.__getattr__('value', filter_ind), axis = 0), np.nanmean(self.n_chg.__getattr__('value', filter_ind), axis = 0), np.nanmean(self.n_chh.__getattr__('value', filter_ind), axis = 0))), columns=["CG","CHG","CHH"] )
        t_req_gene = t_req_gene.set_index( self.cg.accessions )
        t_n_req_gene = t_n_req_gene.set_index( self.cg.accessions )
        return((t_req_gene, t_n_req_gene))

    def get_meths_req_gene(self, req_gene_ix, context_ix = 0, count_thres=10):
        ##  For a given gene return a pandas dataframe methylation and number of cytosine counts
        ## req_genes_str = ["Chr1,10319732,10320402"]
        ## context_ix == 0 ? CG, context_ix == 1 ? CHG, context_ix == 2 ? CHH
        (t_req_gene, t_n_req_gene) = self.get_cg_chg_chh_meths( np.array([req_gene_ix], dtype=int) )
        t_req_gene = t_req_gene.iloc[:, context_ix]
        t_n_req_gene = t_n_req_gene.iloc[:, context_ix]
        t_req_gene[ np.where( t_n_req_gene < count_thres )[0] ] = np.nan
        return(t_req_gene)

    def iterate_meths_req_genes(self, req_genes_str, context_ix=0):
        ## req_genes_str = ["Chr1,10319732,10320402"]
        req_genes_ix = self.get_filter_inds(req_genes_str)
        for egene_ix in req_genes_ix:
            yield( self.get_meths_req_gene(egene_ix, context_ix ) )

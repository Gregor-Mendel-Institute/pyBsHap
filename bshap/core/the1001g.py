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
import itertools
import pybedtools as pybed

from pygenome import genome as g
genome = g.GenomeClass("at_tair10")

log = logging.getLogger(__name__)

def die(msg):
  sys.stderr.write('Error: ' + msg + '\n')
  sys.exit(1)

class WriteHDF51001Table(object):

    def __init__(self, input_file, output_file, chunk_size=1000):
        self.chunk_size = chunk_size
        self.input_file = input_file
        self.output_file = output_file
        # self.genome = g.ArabidopsisGenome(ref_genome)
        #self.num_lines = len(input_file)

    @staticmethod
    def _read_bg_file(bg_file, value_column, header = False, delimiter = "\t"):
        ## Also provide the column which value contains    
        if header:
            e_bg = pd.read_csv(bg_file, header = 0, sep= delimiter,  dtype = {0: str, 1: np.int32, 2: np.int32})
        else:
            e_bg = pd.read_csv(bg_file, header = None, sep= delimiter, dtype = {0: str, 1: np.int32, 2: np.int32})
        e_bg.rename(columns={0:'chr',1:'start',2:'end', value_column - 1:'value' }, inplace=True)
        return(e_bg)

    def parse_accs_ids_input_file(self, split_str = "."):
        ## File names can  generally be
        # 1. mhl.10001.txt
        # 2. 17328_CATTTT_C3N13ACXX_3_20140219B_20140219.snpvcf.genotyper.txt
        input_ids = pd.Series([ os.path.basename(efile) for efile in  self.input_file]).str.split(split_str, expand=True)
        if len(np.unique(input_ids.iloc[:,0])) == len(self.input_file):
            input_ids = np.array(input_ids.iloc[:,0], dtype="str")
        elif len(np.unique(input_ids.iloc[:,1])) == len(self.input_file):
            input_ids = np.array(input_ids.iloc[:,1], dtype="str")
        else:
            nunique = input_ids.apply(pd.Series.nunique)
            cols_to_join = nunique[nunique != 1].index
            input_ids = input_ids[cols_to_join].agg(split_str.join, axis=1)
        return(input_ids.values)

    def write_h5_multiple_files(self, value_column = 4):
        ## self.input_file is an array of files
        log.info("reading %s input files" % len(self.input_file))
        num_lines = len(self.input_file)
        log.info("writing a h5 file for %s bdg files into %s" % (num_lines, self.output_file))
        h5file = h5.File(self.output_file, 'w')
        h5file.create_dataset('files', data=np.array(self.input_file).astype("S"), shape=(num_lines,))
        h5file.create_dataset('accessions', data = self.parse_accs_ids_input_file(), dtype = h5.special_dtype(vlen=bytes) )
        base_bg = self._read_bg_file(self.input_file[0], value_column)
        n_rows = base_bg.shape[0]
        if n_rows < self.chunk_size:
            self.chunk_size = n_rows
        h5file.create_dataset('chunk_size', data=self.chunk_size, shape=(1,),dtype='i8')
        h5file.create_dataset('value', shape = ((n_rows, num_lines)), chunks=((self.chunk_size,num_lines)),dtype='float')
        for ef_ind in range(num_lines - 1):
            e_bg = self._read_bg_file(self.input_file[ef_ind + 1], value_column)
            if n_rows == 0:
                n_rows = e_bg.shape[0]
            else:
                if n_rows != e_bg.shape[0]:
                    die("please provide bedgraph files with exactly same number of lines")
            h5file['value'][:,ef_ind+1] = np.array(e_bg['value'], dtype='float')
            if ef_ind % 50 == 0 and ef_ind > 0:
                log.info("written %s files into h5file" % ef_ind)
        # h5file.create_dataset('chr', shape=(n_rows,), data = np.array([genome.chrs[e] for e in genome.get_chr_ind(base_bg['chr']) ]) )
        h5file.create_dataset('chr', shape=(n_rows,), data = base_bg['chr'].values.astype('S') )
        h5file.create_dataset('start', shape=(n_rows,), data = base_bg['start'].values)
        h5file.create_dataset('end', shape=(n_rows,), data = base_bg['end'].values)
        h5file['value'][:,0] = np.array(base_bg['value'])
        h5file.close()
        log.info("done")

    @staticmethod
    def _read_value_matrix_file(matrix_file):
        v_file = pd.read_csv(matrix_file, header = None, sep= ",")
        return(v_file)

    def write_h5_matrix(self):
        v_file = self._read_value_matrix_file(self.input_file[0])
        num_lines = v_file.shape[1] - 3
        n_rows = v_file.shape[0] - 2
        n_rows_split = v_file.iloc[2:,0].str.split(":", expand =True)
        h5file = h5.File(self.output_file, 'w')
        h5file.create_dataset('files', data=np.array(self.input_file))
        h5file.create_dataset('accessions', data = np.array(v_file.iloc[0,:][3:], dtype="string"), shape= ((num_lines,)), dtype='S20')
        if n_rows < self.chunk_size:
            self.chunk_size = n_rows
        h5file.create_dataset('chunk_size', data=self.chunk_size, shape=(1,),dtype='i8')
        h5file.create_dataset('value', shape = ((n_rows, num_lines)), chunks=((self.chunk_size,num_lines)),dtype='float')
        for ef_ind in range(num_lines):
            e_bg = np.array(v_file.iloc[2:,ef_ind+3], dtype="string")
            e_bg[e_bg == "AA"] = 0
            e_bg[e_bg == "AB"] = 1
            e_bg[e_bg == "BB"] = 2
            h5file['value'][:,ef_ind] = np.array(e_bg, dtype='float')
            if ef_ind % 50 == 0 and ef_ind > 0:
                log.info("written %s files into h5file" % ef_ind)
        h5file.create_dataset('chr', shape=(n_rows,), data = np.array([genome.chrs[e] for e in genome.get_chr_ind(n_rows_split.iloc[:,0]) ]))
        h5file.create_dataset('start', shape=(n_rows,), data = np.array(n_rows_split.iloc[:,1].str.split("-", expand=True).iloc[:,0], dtype='int'))
        h5file.create_dataset('end', shape=(n_rows,), data = np.array(n_rows_split.iloc[:,1].str.split("-", expand=True).iloc[:,1], dtype='int'))
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
        self.accessions = np.array(self.h5file['accessions']).astype('U')

    def __getattr__(self, name, filter_pos_ix=None, return_np=False):
        if name not in ['chr', 'start', 'end', 'value']:
            raise AttributeError("%s is not in the keys for HDF5. Only accepted values are ['chr', 'start', 'end', 'value']" % name)
        if filter_pos_ix is None:
            if return_np:
                if name in ['chr']:
                    return( np.array(pd.Series(self.h5file[str(name)]).str.decode("utf-8") ) )
                return(np.array(self.h5file[str(name)]))
            return(self.h5file[str(name)])
        elif type(filter_pos_ix) is np.ndarray:
            rel_pos_ix = filter_pos_ix - filter_pos_ix[0]
            ret_attr = np.array(self.h5file[str(name)][filter_pos_ix[0]:filter_pos_ix[-1]+1][rel_pos_ix])
        else:
            ret_attr = np.array(self.h5file[str(name)][filter_pos_ix])
        if name in ['chr']:
            ret_attr = np.array( pd.Series(ret_attr).str.decode("utf-8") )
        elif name in ['start', 'end']:
            ret_attr = ret_attr.astype(int)
        else:
            ret_attr = ret_attr.astype(float)
        return(ret_attr)

    def change_accs_ids(self, sra_table):
        ## /projects/cegs/rahul/009.1001methylomes.rawdata/final_dataset_1001meths_rawdata.txt
        sra_table_1001g = pd.read_csv(sra_table, header = None, sep = "\t")
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
        if isinstance(region_bed, str):
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

    def iter_inds_overlap_regions(self, region_bed_df, g = None ):
        assert type(region_bed_df) is pd.core.frame.DataFrame, "please provide a pandas dataframe object"
        self_bed_str = self.get_bed_df(None, return_str = True)
        self_bed = pybed.BedTool.from_dataframe( self.get_bed_df(None) )
        region_bed_str = np.array(region_bed_df.iloc[:,0].map(str) + "," + region_bed_df.iloc[:,1].map(str) + "," +  region_bed_df.iloc[:,2].map(str), dtype="str")
        ## Just taking first three columns for bedtools
        for req_bed_ix in range(region_bed_df.shape[0]):
            t_region_bed = pybed.BedTool.from_dataframe(region_bed_df.iloc[:,[0,1,2]].iloc[[req_bed_ix]] )
            inter_gene_bed = self_bed.intersect(t_region_bed, wa=True)
            if inter_gene_bed.count() == 0:   ## Return if there are no matching lines.
                yield(None)
            inter_gene_bed = inter_gene_bed.to_dataframe()
            inter_gene_bed_str = np.array(inter_gene_bed.iloc[:,0].map(str) + "," + inter_gene_bed.iloc[:,1].map(str) + "," +  inter_gene_bed.iloc[:,2].map(str), dtype="str")
            yield(np.where( np.in1d(self_bed_str, inter_gene_bed_str ) )[0])


    def get_inds_overlap_region_file(self, region_file, just_names=False, araport11_file=None):
        whole_bed = self.get_bed_df(None)
        return(run_bedtools.get_filter_bed_ix(region_file, whole_bed, just_names=just_names, araport11_file=araport11_file))

    def get_inds_matching_region(self, region_bed):
        # returns values given region_bed (a pd dataframe with chr, start and end)
        ## maybe useful for wma hdf5 files
        ## region_bed = ['Chr1',3631, 5899]
        ## or region_bed = "Chr1,3631,5899"
        if isinstance(region_bed, str):
            t_split = region_bed.split(",")
            assert len(t_split) == 3
            region_bed = [ t_split[0], int(t_split[1]), int(t_split[2]) ]
        chr_inds = np.where(self.__getattr__("chr", return_np=True) == region_bed[0])[0]
        chr_start = self.__getattr__("start", chr_inds)
        chr_end = self.__getattr__("end", chr_inds)
        return(chr_inds[np.where( (chr_start == region_bed[1]) & (chr_end == region_bed[2]) )[0]])

    def get_inds_matching_region_file(self, region_file):
        region_df = pd.read_csv(region_file, sep = "\t'")
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
        self._load_files()

    def _load_files(self):
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
        all_strs = np.unique(all_strs, return_index = True)
        req_inds = all_strs[1][np.where( np.in1d(all_strs[0], req_genes_str ) )[0]]
        log.warn("Make sure provided array is sorted based on chromosome position")
        return(np.sort(req_inds))

    def get_average_cg_chg_chh_meths(self, req_gene_ix, count_thres = 10):
        """
        Function to calculate the average CG, CHG and CHH methylation
        input: 
        req_gene_ix: numpy array of bed postion in string, array(['Chr1,17024,18924', 'Chr1,18331,18642', 'Chr1,55676,56576', ...,
       'Chr5,26829819,26831418', 'Chr5,26835012,26835916',
       'Chr5,26891559,26893445'], dtype='<U32')
       count_thres: Filter regions that have smaller than n cytosines called
        """
        ## Given list of indices, function outputs average of methylations in three contexts
        if req_gene_ix is None or len(req_gene_ix) == 0:
            return(None)
        t_meths = self.get_meths_req_gene(req_gene_ix, count_thres)
        t_meths = pd.DataFrame( 
            np.column_stack((
                np.nanmean(t_meths['CG'], axis = 0), 
                np.nanmean(t_meths['CHG'], axis = 0), 
                np.nanmean(t_meths['CHH'], axis = 0) 
            ) ), 
            columns= ['CG', 'CHG', 'CHH'], 
            index= t_meths['CG'].columns
        )
        return( t_meths )

    def differential_methylation(self, accs_ix, filter_pos_ix=None, count_thres=10):
        ## accs_ix = [[0,1], [2,3]]
        assert len(accs_ix) == 2, "provide an array with two classes. eg. [[0,1], [2,3]]"
        meths_all = self.get_meths_req_gene( filter_pos_ix, count_thres = count_thres )
        req_cg = meths_all['CG'].iloc[:,accs_ix[0]].mean(axis = 1) - meths_all['CG'].iloc[:,accs_ix[1]].mean(axis = 1)
        req_chg = meths_all['CHG'].iloc[:,accs_ix[0]].mean(axis = 1) - meths_all['CHG'].iloc[:,accs_ix[1]].mean(axis = 1)
        req_chh = meths_all['CHH'].iloc[:,accs_ix[0]].mean(axis = 1) - meths_all['CHH'].iloc[:,accs_ix[1]].mean(axis = 1)
        return( pd.DataFrame( np.column_stack( ( req_cg, req_chg, req_chh ) ), columns = ['CG', "CHG", 'CHH'] ) )

    def get_meths_req_gene(self, req_gene_ix, count_thres=10):
        ##  For a given gene return a pandas dataframe methylation and number of cytosine counts
        ## req_gene_ix is an np array
        assert type(req_gene_ix) is np.ndarray, "provide an np array for indices"
        t_cg = self.__getattribute__('cg').__getattr__('value', req_gene_ix)
        t_chg = self.__getattribute__('chg').__getattr__('value', req_gene_ix)
        t_chh = self.__getattribute__('chh').__getattr__('value', req_gene_ix)
        t_cg[ self.__getattribute__('n_cg').__getattr__('value', req_gene_ix) < count_thres ] = np.nan
        t_chg[ self.__getattribute__('n_chg').__getattr__('value', req_gene_ix) < count_thres ] = np.nan
        t_chh[ self.__getattribute__('n_chh').__getattr__('value', req_gene_ix) < count_thres ] = np.nan
        req_meths = {'CG': pd.DataFrame(t_cg, columns=self.cg.accessions), 'CHG': pd.DataFrame(t_chg, columns=self.cg.accessions), 'CHH': pd.DataFrame(t_chh, columns=self.cg.accessions)}
        return(req_meths)

    def total_methylations_(self, req_genes_str, outFile=None, count_thres = 5):
        all_meths = self.get_meths_req_gene(self.get_filter_inds(req_genes_str), count_thres = count_thres)
        total_cg_meths = pd.DataFrame( np.column_stack((self.cg.accessions, np.nanmean(np.array(all_meths['CG']), axis = 0) )), columns = ["accessionid", "pheno"] ).dropna()
        total_chg_meths = pd.DataFrame( np.column_stack((self.cg.accessions, np.nanmean(np.array(all_meths['CHG']), axis = 0) )), columns = ["accessionid", "pheno"] ).dropna()
        total_chh_meths = pd.DataFrame( np.column_stack((self.cg.accessions, np.nanmean(np.array(all_meths['CHH']), axis = 0) )), columns = ["accessionid", "pheno"] ).dropna()
        if outFile is not None:
            total_cg_meths.to_csv( outFile + ".CG.txt", index = False )
            total_chg_meths.to_csv( outFile + ".CHG.txt", index = False )
            total_chh_meths.to_csv( outFile + ".CHH.txt", index = False )
        return((total_cg_meths, total_chg_meths, total_chh_meths))

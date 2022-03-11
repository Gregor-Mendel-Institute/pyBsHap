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
from pygenome import genome
import csv
import itertools
import re

log = logging.getLogger(__name__)


def die(msg):
  sys.stderr.write('Error: ' + msg + '\n')
  sys.exit(1)

class writeHDF5MethTable(object):
    """
    Main class function to write hdf5 files from allc, bed files
    """

    def __init__(self, ref_fasta_file, allc_file = None, output_file = None):
        """
        initilize class function
        """
        self.fasta =  genome.GenomeClass(ref_fasta_file)
        if allc_file is not None:
            self.load_allc_file(allc_file)
        if output_file is not None:
            self.write_h5_file( output_file )

    def load_allc_file(self, allc_file, umeth = None):
        log.info("reading the allc file!")
        allc_bed = pd.read_csv(allc_file, sep = "\t", header = None, dtype = {0: "str", 1: np.int64})
        allc_bed.columns = np.array(['chr','pos','strand','mc_class','mc_count','total','methylated'])
        allc_bed_sorted = pd.DataFrame(columns = allc_bed.columns)
        chrpositions = np.zeros(1, dtype="int")
        log.info("sorting the bed file!")
        for ec in self.fasta.chrs:
            t_echr = allc_bed.iloc[np.where(allc_bed['chr'] == ec)[0], :].sort_values(by='pos', ascending=True)
            allc_bed_sorted = allc_bed_sorted.append(t_echr, ignore_index=True)
            chrpositions = np.append(chrpositions, chrpositions[-1] + t_echr.shape[0])
        log.info("done!")
        self.allc_bed = allc_bed_sorted
        self.chrpositions = chrpositions
        if umeth is not None:
            conv_rates = allc_bed.loc[allc_bed['chr'] == umeth, ["mc_count", 'total']].sum(0)
            return( (conv_rates['mc_count'] / conv_rates['total']) )
        return(None)

    def load_bismark_coverage(self, bismark_cov_file, umeth = None, min_total = 1):
        bsbed = pd.read_csv(bismark_cov_file, sep = "\t", header = None, dtype = {0: "str", 1: np.int64})
        bsbed.columns = ['chr', 'pos', 'strand', 'mc_count', 'c_count', 'context', 'mc_class']
        bsbed['total'] = bsbed['mc_count'] + bsbed['c_count']
        self.allc_bed = pd.DataFrame(columns = bsbed.columns)
        self.chrpositions = np.zeros(1, dtype="int")
        for ec in self.fasta.chrs:
            t_echr = bsbed.loc[bsbed['chr'] == ec].sort_values(by='pos', ascending=True)
            t_echr = t_echr[t_echr['total'] >= min_total]
            self.allc_bed = self.allc_bed.append(t_echr, ignore_index=True)
            self.chrpositions = np.append( self.chrpositions, self.chrpositions[-1] + t_echr.shape[0] )
        if umeth is not None:
            conv_rates = bsbed.loc[bsbed['chr'] == umeth, ["mc_count", 'total']].sum(0)
            return( (conv_rates['mc_count'] / conv_rates['total']) )
        return(None) 


    def write_h5_file(self, output_file):
        chunk_size = min(100000, self.allc_bed.shape[0])
        h5file = h5.File(output_file, 'w')
        h5file.create_dataset('chunk_size', data=chunk_size, shape=(1,),dtype='i8')
        h5file.create_dataset('chrpositions', data=self.chrpositions, shape=(len(self.chrpositions),),dtype='int')
        allc_columns = ['chr', 'pos', 'strand', 'mc_class', 'mc_count', 'total']
        allc_dtypes = ['S16', 'int', 'S4', 'S8', 'int', 'int']
        ## Going through all the columns
        for cols, coltype in zip(allc_columns, allc_dtypes):
            h5file.create_dataset(cols, compression="gzip", data=np.array(self.allc_bed[cols],dtype=coltype), shape=(self.allc_bed.shape[0],), chunks = ((chunk_size,)))
        h5file.close()

# def read_allc_files_chrwise(allc_id, allc_path):
#     allcBed = []
#     chrpositions = []
#     tlen = 0
#     sample_chrs = get_chrs(allc_id, allc_path)
#     for c in sample_chrs:
#         allc_file = allc_path + "/allc_" + allc_id + "_" + c + ".tsv"
#         log.info("progress: reading %s!" % allc_file)
#         chrpositions.append(tlen)
#         bsbed = read_allc_pandas_table(allc_file)
#         tlen = tlen + bsbed.shape[0]
#         try:
#             allcBed = pd.concat([allcBed, bsbed])
#         except TypeError:
#             allcBed = bsbed
#     return((allcBed, chrpositions))

def load_hdf5_methylation_file(hdf5_file, ref_fasta="at_tair10", bin_bed=None):
    return HDF5MethTable(hdf5_file, ref_fasta, bin_bed)

# Try to make it as a class, learned from PyGWAS
class HDF5MethTable(object):

    def __init__(self, hdf5_file, ref_fasta = "at_tair10", bin_bed=None):
        self.h5file = h5.File(hdf5_file,'r')
        self.filter_pos_ix = self.get_filter_inds(bin_bed)
        self.chrpositions = np.array(self.h5file['chrpositions'])
        self.genome = genome.GenomeClass(ref_fasta)

    def close_h5file(self):
        if self.h5file is not None:
            self.h5file.close()

    def get_chrinds(self, chrid):
        req_chr_ind = self.genome.get_chr_ind(chrid)
        if req_chr_ind is None:
            return( (None, None) )
        chr_inds = [self.chrpositions[req_chr_ind], self.chrpositions[req_chr_ind + 1]]
        return((req_chr_ind, chr_inds))

    def get_filter_inds(self, bin_bed = None):
        # bin_bed = ['Chr1', 0, 100] or "Chr1,1,100"
        if isinstance(bin_bed, pd.DataFrame):
            t_chr = np.unique(bin_bed.iloc[:,0])
            if t_chr.shape[0] == 1:
                req_chr_ind, chr_inds = self.get_chrinds(t_chr[0])
                refBed = self.get_bed_df(filter_pos_ix = np.arange(chr_inds[0], chr_inds[1]), full_bed=False)
            else:
                refBed = self.get_bed_df(filter_pos_ix = np.arange(self.__getattr__('pos').shape[0]), full_bed=False)
            req_inds_df = run_bedtools.get_intersect_bed_ix(reference_bed=refBed, query_bed=bin_bed, just_names=False)
            if t_chr.shape[0] == 1:
                req_inds_df['ref_ix'] = req_inds_df['ref_ix'] + chr_inds[0]
            return(req_inds_df)
        elif bin_bed is None:
            return(None)
        elif isinstance(bin_bed, str):
            t_split = bin_bed.split(",")
            assert len(t_split) == 3, "please genomic position as 'Chr1,1,100'"
            bin_bed = [ t_split[0], int(t_split[1]), int(t_split[2]) ]
        assert len(bin_bed) == 3, "please genomic position as ['Chr1', 1, 100]"
        req_chr_ind, chr_inds = self.get_chrinds(bin_bed[0])
        if req_chr_ind is None:
            return(np.zeros(0, dtype=int))
        req_inds = chr_inds[0] + np.searchsorted(self.positions[chr_inds[0]:chr_inds[1]], [bin_bed[1],bin_bed[2]], side='right')
        return(np.arange(req_inds[0],req_inds[1]))

    def __getattr__(self, name, filter_pos_ix=None, return_np=False):
        req_names = ['chr', 'position', 'pos', 'positions', 'methylated', 'strand', 'mc_count', 'mc_class', 'total', 'mc_total']
        if name not in req_names:
            raise AttributeError("%s is not in the keys for HDF5. Only accepted values are 'chr', 'positions', 'pos', 'methylated', 'strand', 'mc_class', 'total', 'mc_total'" % name)
        if name in ['total', 'mc_total']:
            name = 'total'
        if name in ['pos', 'position', 'positions']:
            name = 'pos'
        if filter_pos_ix is None:
            if return_np: 
                ret_attr = np.array( self.h5file[str(name)] )
            else:
                return( self.h5file[str(name)] )
        elif type(filter_pos_ix) is np.ndarray:
            if len(filter_pos_ix) == 0:
                ret_attr = np.array(())
            else:
                rel_pos_ix = filter_pos_ix - filter_pos_ix[0]
                ret_attr = np.array(self.h5file[str(name)][filter_pos_ix[0]:filter_pos_ix[-1]+1][rel_pos_ix])
        else:
            ret_attr = np.array(self.h5file[str(name)][filter_pos_ix])
        if name in ['chr', 'mc_class', 'strand']:
            ret_attr = ret_attr.astype('U')
        elif name in ['pos', 'total', 'mc_count', 'methylated']:
            ret_attr = ret_attr.astype(int)
        return(ret_attr)

    def get_bed_df(self, filter_pos_ix, full_bed=False, read_threshold=0):
        req_chrs = self.__getattr__('chr', filter_pos_ix, return_np=True)
        req_pos = self.__getattr__('pos', filter_pos_ix, return_np=True)
        conname = self.__getattr__('mc_class', filter_pos_ix, return_np=True )
        strand = self.__getattr__('strand', filter_pos_ix, return_np=True)
        mc_df = pd.DataFrame(np.column_stack((req_chrs, req_pos, req_pos + 1, conname, strand)), columns=['chr', 'start', 'end', 'mc_class', 'strand'])
        if full_bed:
            # mc_df['permeth'] = self.get_permeths(filter_pos_ix, read_threshold=read_threshold)
            mc_df['mc_count'] = self.__getattr__('mc_count', filter_pos_ix, return_np=True)
            mc_df['mc_total'] = self.__getattr__('mc_total', filter_pos_ix, return_np=True)
        return(mc_df)
        

    def get_req_mc_class_ix(self, req_context, filter_pos_ix):
        if filter_pos_ix is None or len(filter_pos_ix) == 0:
            return(filter_pos_ix)
        if req_context is None:
            if filter_pos_ix is not None:
                return(np.arange(len(filter_pos_ix)))
            else:
                return(None)
        import re
        cor = re.compile(req_context)
        np_vmatch = np.vectorize(lambda x:bool(cor.match(x)))
        return(np.where(np_vmatch(self.__getattr__('mc_class', filter_pos_ix, return_np=True )))[0])

    def get_lowfreq(self, filter_pos_ix=None):
        if filter_pos_ix is None:
            return(np.array(self.h5file['lowfreq']))
        elif type(filter_pos_ix) is np.ndarray:
            rel_pos_ix = filter_pos_ix - filter_pos_ix[0]
            return(self.h5file['lowfreq'][filter_pos_ix[0]:filter_pos_ix[-1]+1][rel_pos_ix])
        else:
            return(self.h5file['lowfreq'][filter_pos_ix])

    def get_permeths(self, filter_pos_ix=None, read_threshold=0):
        # If read_threshold is given
        #   if methylated then we have a value
        #   else the value is 0
        # Else the notation is -1
        # methylated_cs = np.array(self.__getattr__('methylated', filter_pos_ix), dtype="float")
        mc_count = self.__getattr__('mc_count', filter_pos_ix, return_np=True)
        mc_total = self.__getattr__('mc_total', filter_pos_ix, return_np=True)
        # calculated_permeths = np.divide(np.multiply(methylated_cs, mc_count), mc_total)
        calculated_permeths = np.divide(mc_count, mc_total)
        if read_threshold == 0:
            return(calculated_permeths)
        else:
            if filter_pos_ix is None:
                filter_pos_ix = np.arange(len(mc_total))
            permeths = np.repeat(-1.0, len(filter_pos_ix))
            new_filter_pos_ix = np.where(mc_total > read_threshold)[0]
            permeths[new_filter_pos_ix] = calculated_permeths[new_filter_pos_ix]
            return(permeths)

    def iter_chr_windows(self, chrid, window_size):
        req_chr_ind = self.genome.get_chr_ind(chrid)
        return(self.iter_bed_windows([chrid, 1, self.genome.golden_chrlen[req_chr_ind]], window_size))

    def iter_bed_windows(self, required_bed, window_size):
        ## required_bed = ["Chr1", 1, 100]
        filter_pos_ix = self.get_filter_inds(required_bed)
        assert type(filter_pos_ix) is np.ndarray, "please provide a numpy array"
        filter_pos = self.__getattr__('pos', filter_pos_ix, return_np=True)
        ind = 0
        for t in range(required_bed[1], required_bed[2], window_size):
            skipped = True
            result = []
            bin_bed = [int(t), min(required_bed[2], int(t) + window_size - 1)]
            for epos in filter_pos[ind:]:
                if epos >= bin_bed[0]:
                    if epos <= bin_bed[1]:
                        result.append(filter_pos_ix[ind])
                    elif epos > bin_bed[1]:
                        skipped = False
                        yield((bin_bed, result))
                        break
                    ind = ind + 1
            if skipped:
                yield((bin_bed, result))

    def MethylationSummaryStats(self, filter_pos_ix, category = 1, req_context = None, min_depth = 3):
        # meths is already filtered for bin_bed positions
        if filter_pos_ix is None:
            return(np.nan)
        if len(filter_pos_ix) == 0:
            return(np.nan)
        # Filtering for context
        if req_context is not None:
            filter_pos_ix = np.array(filter_pos_ix, dtype=int)[self.get_req_mc_class_ix(req_context, filter_pos_ix)]
            if len(filter_pos_ix) == 0:
                return(np.nan)
        mc_total = self.__getattr__('mc_total', filter_pos_ix, return_np=True)
        mc_count = self.__getattr__('mc_count', filter_pos_ix, return_np=True)
        if np.sum(mc_total) == 0:
            return(np.nan)
        if category == 1:   # weighted mean
            methylated_cs = np.where(mc_total >= min_depth)[0]
            # methylated_cs = np.where(self.__getattr__('methylated', filter_pos_ix, return_np=True) == 1)[0]
            return(np.divide(float(np.sum(mc_count[methylated_cs])), np.sum(mc_total[methylated_cs])))
        elif category == 2: # fraction of methylated positions
            methylated = self.__getattr__('methylated', filter_pos_ix, return_np=True)
            meths_len = np.sum(methylated)
            return(float(meths_len)/len(filter_pos_ix))
        elif category == 3: ## absolute means
            return(np.mean(np.divide(np.array(mc_count, dtype="float"), mc_total)))
        elif category == 4:  ### weighted mean but without making correction for converion.
            return(np.divide(float(np.sum(mc_count)), np.sum(mc_total)))
        else:
            raise(NotImplementedError)

    def _AveMethylation_All_Contexts( self, filter_pos_ix, category = 1 ):
        cg_meth = self.MethylationSummaryStats( filter_pos_ix, category=category, req_context = "CG[ATGC]" )
        chg_meth = self.MethylationSummaryStats( filter_pos_ix, category=category, req_context = "C[ATC]G" )
        chh_meth = self.MethylationSummaryStats( filter_pos_ix, category=category, req_context = "C[ATC][ATC]" )
        return( (cg_meth, chg_meth, chh_meth) )

    def _TotalCounts_All_Contexts(self, filter_pos_ix):
        t_count_cg = self.get_req_mc_class_ix("CG[ATGC]", filter_pos_ix)
        t_count_chg = self.get_req_mc_class_ix("C[ATC]G", filter_pos_ix)
        t_count_chh = self.get_req_mc_class_ix("C[ATC][ATC]", filter_pos_ix)
        return( (t_count_cg, t_count_chg, t_count_chh) )

    def generate_meth_average_in_windows(self, out_file, window_size, category=1):
        ## Methylation category here by default is weighted average
        # look the function in HDF51001gTable MethylationSummaryStats).
        outmeths_cg_avg = open(out_file + ".CG.bg", 'w')
        outmeths_chg_avg = open(out_file + ".CHG.bg", 'w')
        outmeths_chh_avg = open(out_file + ".CHH.bg", 'w')
        for echr, echrlen in zip(self.genome.chrs, self.genome.golden_chrlen):
            self_windows = self.iter_chr_windows(echr, window_size)
            count = 0
            log.info("analyzing %s" % echr)
            for ewin in self_windows:
                t_count = self._TotalCounts_All_Contexts( ewin[1] )
                req_meth_avg = self._AveMethylation_All_Contexts(ewin[1], category)
                outmeths_cg_avg.write("%s\t%s\t%s\t%s\t%s\n" % (echr, ewin[0][0], ewin[0][1], req_meth_avg[0], len(t_count[0])))
                outmeths_chg_avg.write("%s\t%s\t%s\t%s\t%s\n" % (echr, ewin[0][0], ewin[0][1], req_meth_avg[1], len(t_count[1])))
                outmeths_chh_avg.write("%s\t%s\t%s\t%s\t%s\n" % (echr, ewin[0][0], ewin[0][1], req_meth_avg[2], len(t_count[2])))
                count = count + 1
                if count % 5000 == 0:
                    log.info("progress: analysed %s windows" % count)
        outmeths_cg_avg.close()
        outmeths_chg_avg.close()
        outmeths_chh_avg.close()
        log.info("done!")

    def generate_meth_average_required_bed(self, required_bed, out_file = None, sort_bed = True, category=1):
        indexed_bed = False
        if type(required_bed) is str:
            if op.isfile(required_bed):
                req_regions = pd.read_csv(required_bed, sep = "\t", header=None)
                if len(req_regions.columns) >= 4:
                    req_regions = req_regions.set_index(3)
                    indexed_bed = True
            else:
                ### required_bed = "Chr1,1,1000"
                req_regions = pd.DataFrame(required_bed.split(",")).T
        elif type(required_bed) is pd.DataFrame:
            req_regions = required_bed
        ## Sort the given bed file
        if sort_bed:
            req_regions = run_bedtools.sort_bed_df(req_regions)
        assert len(req_regions.columns) >= 3, "bed file should have minimum of 3 columns. eg., Chr, Start, End"
        if out_file is not None:
            outmeths_cg_avg = open(out_file + ".CG.bg", 'w')
            outmeths_chg_avg = open(out_file + ".CHG.bg", 'w')
            outmeths_chh_avg = open(out_file + ".CHH.bg", 'w')
        output_meths = pd.DataFrame( index = req_regions.index, columns = ['cg', 'chg', 'chh'] )
        calc_for_each = True
        if req_regions.shape[0] > 10:
            mat_positions = self.get_filter_inds( req_regions )
            calc_for_each = False
        for er in req_regions.iterrows():
            if calc_for_each:
                t_filter_pos_ix = self.get_filter_inds( [er[1].iloc[0], int(er[1].iloc[1]), int(er[1].iloc[2]) ] )
            else:
                t_filter_pos_ix = mat_positions.loc[mat_positions['query_ix'] == np.where(req_regions.index.values == er[0])[0][0],'ref_ix'].values
            count = 0
            t_count = self._TotalCounts_All_Contexts( t_filter_pos_ix )
            req_meth_avg = self._AveMethylation_All_Contexts(t_filter_pos_ix, category)
            if out_file is not None:
                outmeths_cg_avg.write("%s\t%s\t%s\t%s\t%s\n" % (er[1].iloc[0], int(er[1].iloc[1]), int(er[1].iloc[2]), req_meth_avg[0], len(t_count[0])))
                outmeths_chg_avg.write("%s\t%s\t%s\t%s\t%s\n" % (er[1].iloc[0], int(er[1].iloc[1]), int(er[1].iloc[2]), req_meth_avg[1], len(t_count[1])))
                outmeths_chh_avg.write("%s\t%s\t%s\t%s\t%s\n" % (er[1].iloc[0], int(er[1].iloc[1]), int(er[1].iloc[2]), req_meth_avg[2], len(t_count[2])))
            output_meths.loc[er[0], 'cg'] = req_meth_avg[0]
            output_meths.loc[er[0], 'chg'] = req_meth_avg[1]
            output_meths.loc[er[0], 'chh'] = req_meth_avg[2]
            count = count + 1
            if count % 100 == 0:
                log.info("progress: analysed %s regions" % count)
        if out_file is not None:
            outmeths_cg_avg.close()
            outmeths_chg_avg.close()
            outmeths_chh_avg.close()
            if indexed_bed:
                output_meths.to_csv( out_file + ".csv" )
        log.info("done!")
        return(output_meths)

    def calculate_gbm_exon_only(self, input_gff_db, gene_id, out_file = None):
        """
        Function to calculate gene-body methylation for a given gene ID
        input:
            input_gff_db:  DB file from gffutils
            gene_id: pandas series for all the IDs
        """
        assert type(gene_id) == pd.Series, "please provide a pandas series for the annotation IDs"
        t_gene_exons = pd.DataFrame( columns=['chr','start', 'end', 'geneid'] )
        gene_start_end = pd.DataFrame(columns=['chr','start', 'end'], index = gene_id )
        for ef_gene in gene_id:
            gene_start_end.loc[ef_gene,:] = [input_gff_db[ef_gene].chrom, input_gff_db[ef_gene].start, input_gff_db[ef_gene].end ]
            for i in input_gff_db.children(input_gff_db[ef_gene], featuretype='exon', order_by='start'):
                t_gene_exons.loc[ i.id, 'geneid' ] = ef_gene
                t_gene_exons.loc[ i.id, 'chr' ] = i.chrom
                t_gene_exons.loc[ i.id, 'start' ] = i.start
                t_gene_exons.loc[ i.id, 'end' ] = i.end
        t_gene_exons['start'] = t_gene_exons['start'].astype(int)
        t_gene_exons['end'] = t_gene_exons['end'].astype(int)
        all_gene_meths = pd.DataFrame( columns=['cg','chg', 'chh'], index = gene_id )
        cg_gene_meths = pd.concat([gene_start_end, pd.DataFrame( columns=['cg','ncg'], index = gene_id )], axis = 1)
        chg_gene_meths = pd.concat([gene_start_end, pd.DataFrame( columns=['chg','nchg'], index = gene_id )], axis = 1)
        chh_gene_meths = pd.concat([gene_start_end, pd.DataFrame( columns=['chh','nchh'], index = gene_id )], axis = 1)
        if len(gene_id) < 100:
            t_filter_pos_ix = np.zeros(0, dtype = int)
            for ef_gene in gene_id:
                ef_gene_exon = t_gene_exons.iloc[ np.where(t_gene_exons['geneid'] == ef_gene)[0],: ]
                t_gene_cds_pos_ix = np.zeros(0, dtype = int)
                for ef_exon in ef_gene_exon.iterrows():
                    t_gene_cds_pos_ix = np.append(t_gene_cds_pos_ix, self.get_filter_inds( ef_exon[1].values[0:3] ))
                t_gene_cds_pos_ix = np.unique(t_gene_cds_pos_ix)
                ef_gene_meths = list(self._AveMethylation_All_Contexts(t_gene_cds_pos_ix, 1))
                ef_gene_meths_counts = list(map(len, self._TotalCounts_All_Contexts( t_gene_cds_pos_ix )))
                all_gene_meths.loc[ef_gene,:] = ef_gene_meths
                cg_gene_meths.loc[ ef_gene, ['cg', 'ncg'] ] = [ef_gene_meths[0], ef_gene_meths_counts[0]]
                chg_gene_meths.loc[ ef_gene, ['chg', 'nchg'] ] = [ef_gene_meths[1], ef_gene_meths_counts[1]]
                chh_gene_meths.loc[ ef_gene, ['chh', 'nchh'] ] = [ef_gene_meths[2], ef_gene_meths_counts[2]]
        else: 
            all_genes_cds_pos_ix = self.get_filter_inds( t_gene_exons )
            for ef_gene in gene_id:
                t_gene_cds_pos_ix = np.unique( 
                    all_genes_cds_pos_ix.loc[
                        np.isin(
                            all_genes_cds_pos_ix['query_ix'], 
                            np.where( t_gene_exons['geneid'] == ef_gene  )[0]
                        ), 
                        'ref_ix'
                    ].values)
                ef_gene_meths = list(self._AveMethylation_All_Contexts(t_gene_cds_pos_ix, 1))
                ef_gene_meths_counts = list(map(len, self._TotalCounts_All_Contexts( t_gene_cds_pos_ix )))
                all_gene_meths.loc[ef_gene,:] = ef_gene_meths
                cg_gene_meths.loc[ ef_gene, ['cg', 'ncg'] ] = [ef_gene_meths[0], ef_gene_meths_counts[0]]
                chg_gene_meths.loc[ ef_gene, ['chg', 'nchg'] ] = [ef_gene_meths[1], ef_gene_meths_counts[1]]
                chh_gene_meths.loc[ ef_gene, ['chh', 'nchh'] ] = [ef_gene_meths[2], ef_gene_meths_counts[2]]
        if out_file is not None:
            cg_gene_meths.to_csv( out_file + ".CG.bg", sep = "\t", index = False, header = None )
            chg_gene_meths.to_csv( out_file + ".CHG.bg", sep = "\t", index = False, header = None )
            chh_gene_meths.to_csv( out_file + ".CHH.bg", sep = "\t", index = False, header = None )
        return(all_gene_meths)

    # def calculate_average_deviations_from_fixed(self, bed_file = None):



def potatoskin_methylation_averages(args):
    # bin_bed = Chr1,1,100
    if args['allc_path'] == "bed":
        args['outFile'] =  args['outFile'] + '.bedGraph'
        run_bedtools.get_genomewide_methylation_average(args)
        log.info("finished!")
        return(0)
    meths = load_hdf5_methylation_file(args['inFile'], args['ref_genome']) 
    if args['required_region'] == '0,0,0':
        meths.generate_meth_average_in_windows(args['outFile'], args['window_size'], category=args['category'])
        return(0)
    meths.generate_meth_average_required_bed(args['required_region'], args['outFile'], args['category'])
    log.info('finished!')
    return(0)


def potatoskin_calculate_gbm_exon_only(args):
    import gffutils
    meths = load_hdf5_methylation_file(args['inFile'], args['ref_genome']) 
    ## Generate a gffutils DB from the a GFF annotation file (check documentation for gffutils)
    ## Assuming 4th column is the gene ID from a bed file
    gff_db = gffutils.FeatureDB(args['gffutils_db'], keep_order=True)
    gene_bed = pd.read_csv( args['gene_bed'], header = None, sep = "\t" )
    t = meths.calculate_gbm_exon_only(gff_db, pd.Series(gene_bed.iloc[:,3].values), args['outFile'])
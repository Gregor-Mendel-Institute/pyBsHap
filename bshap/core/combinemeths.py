# Main module for bisulfite analysis
# Summary statistics
from doctest import OutputChecker
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
import re
import itertools

log = logging.getLogger(__name__)

from . import meth5py
from pygenome import genome as g

def np_get_fraction(x, y, y_min = 0):
    with np.errstate(divide='ignore', invalid='ignore'):
        p = np.divide( x, y )
        if type(y) is np.ndarray:
            p[np.where(y <= y_min)] = np.nan
        else:
            if y <= y_min:
                p = np.nan            
    return(p)

def getOptimumChunks(num_cols, bits_per_dtype = 4, max_mb_io_speed = 10):
    ## IO is the basic estimate.
    ## Maximum chunk size is 4Gb. 
    ## here I am using chunk size to be approximately 1Gb
    # Get bits per dtype as np.float32().itemsize
    perrow_mb = num_cols * bits_per_dtype #* 1024 * 1024
    num_rows = int((max_mb_io_speed * 1024 * 1024) / perrow_mb)
    return( (int(num_rows / 1000) + 1) * 1000 )

class CombinedMethsTable(object):
    ## class for a combined meths tables

    def __init__(self, file_paths, file_ids = None, genome = "at_tair10", load_all= True):
        if load_all:
            self.meths_list = self.load_meths_files(file_paths)
        self._file_path = file_paths
        self.num_lines = len(file_paths)
        self.genome = g.GenomeClass(genome)
        if file_ids is not None:
            self.file_ids = file_ids
        else:
            self.file_ids = np.array([op.basename( ef ) for ef in file_paths])

    def load_meths_files(self, file_paths):
        log.info("reading input files")
        meths_list = []
        for e_path in file_paths:
            meths_list.append(meth5py.load_hdf5_methylation_file(e_path))
        log.info("done! input %s files..!" % len(meths_list))
        return(meths_list)

    def get_cytosines_df_in_bed(self, bed_str, mc_context = "CGN"):
        bed_str = bed_str.split(",")
        bed_str = [bed_str[0], int(bed_str[1]), int(bed_str[2])]
        mc_permeths_df = self.meths_list[0].get_bed_df( self.meths_list[0].get_filter_inds( bed_str ), full_bed = True, read_threshold = 0 )
        if mc_permeths_df.shape[0] == 0:
            mc_permeths_df = self.meths_list[1].get_bed_df( self.meths_list[1].get_filter_inds( bed_str ), full_bed = True, read_threshold = 0 )
        mc_permeths_df = mc_permeths_df.loc[:,['chr', 'start', 'mc_class'] ]
        # mc_permeths_df.columns = ['chr', 'start', 'mc_class', file_ids[0]]
        for m_ix in range(self.num_lines):
            t_meths_df = self.meths_list[m_ix].get_bed_df( self.meths_list[m_ix].get_filter_inds( bed_str ), full_bed = True, read_threshold = 0 )
            t_meths_df = t_meths_df.loc[:,['chr', 'start', 'mc_class', 'mc_count', 'mc_total']]
            # t_meths_df = t_meths_df.loc[t_meths_df['permeth'].astype(float) >= 0, ['chr', 'start', 'mc_class', 'permeth']]
            t_meths_df.columns = ['chr', 'start', 'mc_class', 'mc_count.' + self.file_ids[m_ix], 'mc_total.' + self.file_ids[m_ix]]
            mc_permeths_df = mc_permeths_df.merge( t_meths_df, on = ['chr', 'start', 'mc_class'], how = "outer")
        return( mc_permeths_df )

    def derive_most_common_positions_echr(self, bed_str, num_lines_with_info = 2):
        # returns common positions for a single chromosome
        #from tempfile import mkdtemp
        assert type(bed_str) is str, "provide a str object"
        bed_str = bed_str.split(',')
        bed_str = [str(bed_str[0]), int(bed_str[1]), int(bed_str[2])]
        common_positions = np.arange(bed_str[1], bed_str[2] + 1, dtype="int")
        common_positions_scores = np.zeros( len(common_positions), dtype="int16" )
        t_echr_pos_ix = np.repeat(-1, len(common_positions) * self.num_lines).reshape((len(common_positions), self.num_lines))
        for m_ix in range(self.num_lines):
            m = self.meths_list[m_ix]
            e_chr_pos_ix = m.get_filter_inds(bed_str)
            e_req_pos = m.__getattr__("positions", e_chr_pos_ix, return_np=True)
            common_positions_scores[e_req_pos - bed_str[1]] += 1
            t_echr_pos_ix[e_req_pos - bed_str[1], m_ix] = e_chr_pos_ix
        ## Now calculate which positions are needed
        req_pos_ix = np.where( common_positions_scores >= num_lines_with_info )[0]
        common_positions = common_positions[req_pos_ix]
        return((common_positions, t_echr_pos_ix[req_pos_ix, :] ))

    def derive_most_common_positions(self, output_file, min_mc_total = 3, num_lines_with_info = 2):
        ## derive most common positions for all the chromosomes
        if os.path.exists(output_file):
            return( h5.File(output_file, 'r' ) )
        common_chrs = np.zeros(0, dtype = "S" )
        common_positions = np.zeros(0, dtype = "int" )
        filter_pos_ixs = np.zeros(shape = (0, self.num_lines), dtype = "int")
        for e_chr in self.genome.chrs:
            log.info("reading in through chromosome %s" % e_chr)
            e_bed_str = str(e_chr) + ",1," + str(self.genome.golden_chrlen[self.genome.get_chr_ind(e_chr)])
            e_common_pos, e_filter_pos_ix = self.derive_most_common_positions_echr(e_bed_str, num_lines_with_info)
            common_chrs = np.append(common_chrs, np.repeat( e_chr, e_common_pos.shape[0] ) )
            common_positions = np.append(common_positions, e_common_pos )
            filter_pos_ixs = np.append(filter_pos_ixs, e_filter_pos_ix, axis = 0)
        log.info("done! total number of positions: %s" % str(len(common_positions)))
        self.write_methylated_positions(common_chrs, common_positions, filter_pos_ixs, output_file, min_mc_total = min_mc_total)
        return( h5.File(output_file, 'r' ) )

    def write_methylated_positions(self, common_chrs, common_positions, filter_pos_ixs, output_file, min_mc_total = 3):
        num_positions = common_positions.shape[0]
        chunk_size = min(100000, num_positions )
        # chunk_size = max(100000, min( getOptimumChunks(self.num_lines, 4), num_positions ) )
        log.info("writing the data into h5 file")
        outh5file = h5.File(output_file, 'w')
        data_mc_class = np.repeat( 'nan', num_positions )
        outh5file.create_dataset('chunk_size', data=chunk_size, shape=(1,),dtype='i8')
        outh5file.create_dataset('num_positions', data=num_positions, shape=(1,),dtype='i4')
        outh5file.create_dataset('file_ids', data = np.array(self.file_ids, dtype = "S"))
        outh5file.create_dataset('chr', data = np.array(common_chrs, dtype = "S"), compression = "lzf")
        outh5file.create_dataset('start', dtype = int, data = common_positions, compression = "lzf")
        outh5file.create_dataset('end', dtype = int, data = common_positions + 1, compression = "lzf")
        outh5file.create_dataset('filter_pos_ix', dtype = np.int64, data = filter_pos_ixs, chunks=((chunk_size, self.num_lines)), compression = "lzf")   ### does automatic chunking
        mcs_count = outh5file.create_dataset('mc_count', fillvalue = 0, shape = (num_positions, self.num_lines), dtype = np.int32, chunks=((chunk_size, self.num_lines)), compression = "lzf")
        mcs_total = outh5file.create_dataset('mc_total', fillvalue = 0, shape = (num_positions, self.num_lines), dtype = np.int32, chunks=((chunk_size, self.num_lines)), compression = "lzf")
        mcs_permeths = outh5file.create_dataset('permeths', fillvalue = np.nan, shape = (num_positions, self.num_lines), dtype = np.float32, chunks=((chunk_size, self.num_lines)), compression = "lzf")
        for ef_pos_ix in range(0, num_positions, chunk_size):
            end_point_ix = min(ef_pos_ix+chunk_size, num_positions)
            ef_filter_pos_ixs = filter_pos_ixs[ef_pos_ix:end_point_ix,: ]
            t_ebed_mcs_count = np.zeros( shape = ef_filter_pos_ixs.shape, dtype = int )
            t_ebed_mcs_total = np.zeros( shape = ef_filter_pos_ixs.shape, dtype = int )
            t_ebed_mcs_permeths = np.array(np.tile( np.nan, ef_filter_pos_ixs.shape), dtype = float )
            for ef_mfile in range( self.num_lines ):
                ef_mfile_no_nan = np.where( ef_filter_pos_ixs[:,ef_mfile] != -1 )[0]
                if len(ef_mfile_no_nan) > 0:
                    ef_ref_pos_ix = np.arange( ef_pos_ix, end_point_ix )[ef_mfile_no_nan]
                    ef_mfile_mc_df = self.meths_list[ef_mfile].get_bed_df(ef_filter_pos_ixs[ef_mfile_no_nan,ef_mfile], full_bed = True, read_threshold = 0 )
                    t_ebed_mcs_count[ef_mfile_no_nan, ef_mfile] = ef_mfile_mc_df['mc_count'].values
                    t_ebed_mcs_total[ef_mfile_no_nan, ef_mfile] = ef_mfile_mc_df['mc_total'].values
                    t_ebed_mcs_permeths[ef_mfile_no_nan, ef_mfile] = np_get_fraction(ef_mfile_mc_df['mc_count'].values, ef_mfile_mc_df['mc_total'].values, y_min = min_mc_total  )
                    if np.sum(data_mc_class[ef_ref_pos_ix] == "nan") > 0:
                        data_mc_class[ef_ref_pos_ix] = ef_mfile_mc_df['mc_class'].values
            mcs_count[ef_pos_ix:end_point_ix,:] = t_ebed_mcs_count
            mcs_total[ef_pos_ix:end_point_ix,:] = t_ebed_mcs_total
            mcs_permeths[ef_pos_ix:end_point_ix,:] = t_ebed_mcs_permeths
        outh5file.create_dataset('mc_class', data = np.array(data_mc_class, dtype = "S"))
        outh5file.close()
        log.info("done")

    def derive_methylated_identical_pos_ix_echr(self, chrid, read_depth, pos_in_atleast, max_read_depth = 80):
        ## Caution: would probably need much memory if you give many files.
        # meth_list is a list of meths meth5py object read
        # returns common positions.
        chrid_ind = self.genome.get_chr_ind(chrid)
        num_positions = self.genome.golden_chrlen[chrid_ind]
        common_positions_scores = np.zeros( num_positions, dtype="int" )
        t_echr_pos_ix = np.repeat(-1, num_positions * self.num_lines).reshape((num_positions, self.num_lines))
        for m_ix in range(self.num_lines):
            m = self.meths_list[m_ix]
            e_chrinds = m.get_chrinds(chrid)
            e_chr_pos = m.positions[e_chrinds[1][0]:e_chrinds[1][1]]  ## these are the indices here
            e_chr_depth = m.__getattr__('total', np.arange(e_chrinds[1][0], e_chrinds[1][1]) )
            e_chr_depth[e_chr_depth < read_depth] = 0
            e_chr_depth[e_chr_depth > max_read_depth] = 0
            e_chr_depth[e_chr_depth >= read_depth] = 1
            common_positions_scores[e_chr_pos - 1] = np.add( common_positions_scores[e_chr_pos - 1], m.__getattr__('methylated', np.arange(e_chrinds[1][0], e_chrinds[1][1]) ) * e_chr_depth )
            t_echr_pos_ix[e_chr_pos - 1, m_ix] = np.arange(len(e_chr_pos)) + e_chrinds[1][0]
        req_pos_ix = np.where( common_positions_scores >= pos_in_atleast )[0]
        for m_ix in range(self.num_lines):
            m = self.meths_list[m_ix]
            m_filter_inds = t_echr_pos_ix[req_pos_ix, m_ix]
            if m.filter_pos_ix is None:  ## This is to mainly append when checking for all the chromosomes.
                m.filter_pos_ix = m_filter_inds
            else:
                m.filter_pos_ix = np.append(m.filter_pos_ix, m_filter_inds)
        ### meths_list now has the filter_pos_ix as a key

    def derive_methylated_identical_pos_ix(self, read_depth=5, pos_in_atleast=1, max_read_depth = 80):
        ## Here you get all the common positions going in a loop for each chromosome.
        # It might be some time consuming.
        for e_chr in self.genome.chrs:
            # the function below is appending the list m.filter_pos_ix
            # So, make sure you have some place in the RAM.
            log.info("reading in through chromosome %s" % e_chr)
            self.derive_methylated_identical_pos_ix_echr(e_chr, read_depth=read_depth, pos_in_atleast=pos_in_atleast, max_read_depth=max_read_depth)
        log.info("done!")


class EpiMutations(CombinedMethsTable):
    """
    Class object for estimating epimutation rate

    """
    def __init__(self, file_paths, combined_meths_hdf5_file, file_ids = None, genome = "at_tair10"):
        super().__init__( file_paths, file_ids = file_ids, genome = genome, load_all = False )
        self.f_mcs = h5.File(combined_meths_hdf5_file, 'r')
        self.file_ids = np.array(self.f_mcs['file_ids']).astype("U")

    def search_file_ids(self, search_keys, **kwargs):
        """
        Function to give index of the file ID which contains a specific string
        input:
            provide a pandas series with all the search string. index corresponds to the key
        """
        search_indices = {}
        assert type(search_keys) is pd.Series, "provide a series with search strings and index as key"
        for ef_dict in search_keys.items():
            ef_acc_ix = np.where(pd.Series(self.file_ids).str.contains( ef_dict[1], kwargs))[0]
            if ef_dict[0] in search_indices.keys():
                search_indices[ef_dict[0]] = np.sort(np.append(search_indices[ef_dict[0]], ef_acc_ix))
            else:
                search_indices[ef_dict[0]] = ef_acc_ix
        return(search_indices)

    def __getattr__(self, name, filter_pos_ix=None, return_np=False):
        req_names = ['chr', 'start', 'end', 'count', 'mc_count', 'total', 'mc_total', 'permeths', 'mc_class']
        if name not in req_names:
            raise AttributeError("%s is not in the keys for HDF5. Only accepted values are %s" % (name, req_names))
        if name in ['total', 'mc_total']:
            name = 'mc_total'
        if name in ['count', 'mc_count']:
            name = 'mc_count'
        if filter_pos_ix is None:
            if return_np: 
                ret_attr = np.array( self.f_mcs[str(name)] )
            else:
                return( self.f_mcs[str(name)] )
        elif type(filter_pos_ix) is np.ndarray:
            if len(filter_pos_ix) == 0:
                ret_attr = np.array(())
            else:
                rel_pos_ix = filter_pos_ix - filter_pos_ix[0]
                ret_attr = np.array(self.f_mcs[str(name)][filter_pos_ix[0]:filter_pos_ix[-1]+1][rel_pos_ix])
        else:
            ret_attr = np.array(self.f_mcs[str(name)][filter_pos_ix])
        if name in ['chr', 'mc_class']:
            ret_attr = ret_attr.astype('U')
        # elif name in ['start', 'total', 'mc_count', 'methylated']:
            # ret_attr = ret_attr.astype(int)
        return(ret_attr)

    def get_bed_pos_ix(self, filter_pos_ix, updown = None):
        
        bed_df = pd.DataFrame({
            "chr": self.__getattr__("chr", filter_pos_ix ),
            "start": self.__getattr__("start", filter_pos_ix ),
            "end": self.__getattr__("end", filter_pos_ix )
        }, index = filter_pos_ix)
        if updown is not None:
            bed_df['start'] = bed_df['start'].values - updown
            bed_df['end'] = bed_df['end'].values + updown
        return(bed_df)
    
    def get_req_pos_bed_str(self, bed_str, req_mc_class = "CG[ATCG]", exon_bed_df = None):
        assert type(bed_str) is str, "provide a str object"
        bed_str_split = bed_str.split(',')
        bed_str_split = [str(bed_str_split[0]), int(bed_str_split[1]), int(bed_str_split[2])]
        ## Filter positions that are in a required bed region
        chr_ind_start = np.searchsorted( np.array(self.f_mcs['chr']), np.bytes_(bed_str_split[0]), "left" )
        chr_ind_end = np.searchsorted( np.array(self.f_mcs['chr']), np.bytes_(bed_str_split[0]), "right" )
        req_pos_start = np.arange(chr_ind_start, chr_ind_end)[np.searchsorted(self.f_mcs['start'][chr_ind_start:chr_ind_end], bed_str_split[1]  )]
        req_pos_end = np.arange(chr_ind_start, chr_ind_end)[np.searchsorted(self.f_mcs['start'][chr_ind_start:chr_ind_end], bed_str_split[2]  )]
        ## Filter positions which are in CG context
        filter_pos_ix = np.arange( req_pos_start, req_pos_end )[ np.where(pd.Series(np.array(self.f_mcs['mc_class'][req_pos_start:req_pos_end]).astype("U")).str.contains(req_mc_class, regex = True))[0] ]
        ## Filter CG positions present in gene exons
        if exon_bed_df is not None:
            req_mc_bed = pd.DataFrame(
                {
                    "chr": bed_str_split[0],
                    "start": self.__getattr__('start', filter_pos_ix)
                }
            )
            req_exon_inds = run_bedtools.intersect_positions_bed(reference_bed=exon_bed_df, query_bed=req_mc_bed)
            filter_pos_ix = filter_pos_ix[req_exon_inds]

        return(filter_pos_ix)

    def get_req_pos_ix_genome(self, req_bed_df_dict = None, return_contexts = True, cache_file = None):
        mc_data = {}
        if return_contexts:
            context_search = { 'mcs_mcg_inds': b'CG[ATGC]', 'mcs_mchg_inds': b'C[ATC]G','mcs_mchh_inds': b'C[ATC][ATC]' }
            for ef_context in context_search.keys():
                mc_data[ef_context] = cache_pd_variable_to_file(var_key = ef_context, cache_file = cache_file )
                if mc_data[ef_context] is None:
                    mc_data[ef_context] = pd.Series( np.where( pd.Series(self.f_mcs['mc_class']).apply( lambda x: re.match( context_search[ef_context], x ) is not None ) )[0] )

        if req_bed_df_dict is None:
            return(mc_data)

        query_bed = pd.DataFrame({"chr": np.array(self.f_mcs['chr']).astype('U'), "start": np.array(self.f_mcs['start']) })
        for ef_dict in req_bed_df_dict.keys():
            mc_data['mcs_' + ef_dict + "_inds"] = cache_pd_variable_to_file(var_key = 'mcs_' + ef_dict + "_inds", cache_file = cache_file )
            if mc_data['mcs_' + ef_dict + "_inds"] is None:
                mc_data['mcs_' + ef_dict + "_inds"] = pd.Series( run_bedtools.intersect_positions_bed(reference_bed=req_bed_df_dict[ef_dict], query_bed=query_bed) )
        return(mc_data)


    def calculate_per_meths_per_population(self, sub_populations, filter_cg_pos_ix, max_meth_for_gain = 0.2, min_meth_for_loss = 0.8 , mc_total_min = 3, y_min = 20):
        assert type(sub_populations) is dict, "provide a dictionary with subpoulation ID and index"

        mc_count = self.__getattr__('mc_count', filter_cg_pos_ix ) 
        mc_total = self.__getattr__('mc_total', filter_cg_pos_ix ) 
        # mc_permeth = self.__getattr__('permeths', filter_cg_pos_ix ) 

        epi_out = calculate_deviations_per_populations(
            mc_count, 
            mc_total, 
            sub_populations, 
            max_meth_for_gain = max_meth_for_gain, 
            min_meth_for_loss = min_meth_for_loss, 
            mc_total_min = mc_total_min, 
            prop_y_min = y_min
        )
        epi_out['permeths_subpop'].index = filter_cg_pos_ix
        epi_out['deviations'].index = self.file_ids[epi_out['deviations'].index.values]
        return( epi_out )


def calculate_deviations_per_populations(mc_count, mc_total, sub_populations, max_meth_for_gain = 0.2, min_meth_for_loss = 0.8, mc_total_min = 3, prop_y_min = 20 ):
    """
    Function to estimate the deviations for a given sites. Here is the algorithm
    Input: 
        - Provide 2d numpy array for number of cytosines called for each cytosine
        - provide sub populations as a dictionary with indices (column index in mc_count and mc_total arrays)
    Analysis
        1. Calculate methylation average by pooling inds reads on each cytosine for every subpopulations 
        2. Determine the sites which inherit either a 0 or 1 based on the given thresholds on pooled sub population methylation levels
        3. Calculate gains and loss for each individual at these sites

    Output:
        - Methylation average for each of the given cytosine in each subpopulation
        - Deviations for each individual.
    """
    num_snps = mc_count.shape[0]
    assert mc_count.shape == mc_total.shape, "provide two arrays with same shape"
    mc_permeth = np_get_fraction(mc_count, mc_total, y_min = mc_total_min) 
    epi_out = {}
    epi_out['permeths_subpop'] = pd.DataFrame(index = np.arange( num_snps ) )

    epi_out['deviations'] = pd.DataFrame( columns=['subpop', 'deviation_0', 'deviation_1', 'mc_total_0', 'mc_total_1', 'site_deviation_0', 'site_deviation_1', 'site_total_0', 'site_total_1'] )
    for ef_pop in sub_populations.items():
        epi_out['permeths_subpop'][ef_pop[0]] = np_get_fraction(mc_count[:,ef_pop[1]].sum(1), mc_total[:,ef_pop[1]].sum(1), y_min = prop_y_min)
        t_denovo_ix = np.where(epi_out['permeths_subpop'][ef_pop[0]] <= max_meth_for_gain)[0]
        t_demeth_ix = np.where(epi_out['permeths_subpop'][ef_pop[0]] >= min_meth_for_loss)[0]

        epi_out['permeths_subpop']['inherit_0_' + ef_pop[0]] = 0
        epi_out['permeths_subpop'].loc[t_denovo_ix, 'inherit_0_' + ef_pop[0]]= 1
        epi_out['permeths_subpop']['inherit_1_' + ef_pop[0]] = 0
        epi_out['permeths_subpop'].loc[t_demeth_ix, 'inherit_1_' + ef_pop[0]] = 1

        ## Calculate number of cytosines that increase its methylation
        t_mc_permeths = np.array(mc_permeth[:,ef_pop[1]]).copy()
        # t_mc_permeths[np.greater_equal(t_mc_permeths, 0.5, where = ~np.isnan(t_mc_permeths))] = 1
        # t_mc_permeths[np.less(t_mc_permeths, 0.5, where = ~np.isnan(t_mc_permeths))] = 0
        with np.errstate(invalid='ignore'):
            t_mc_permeths[t_mc_permeths >= 0.5] = 1
            t_mc_permeths[t_mc_permeths < 0.5] = 0
        t_mc_total_0 = (~np.isnan(t_mc_permeths[t_denovo_ix])).sum(0)
        t_mc_total_1 = (~np.isnan(t_mc_permeths[t_demeth_ix])).sum(0)
        t_deviation_0 = np_get_fraction(np.nansum(t_mc_permeths[t_denovo_ix], axis = 0), t_mc_total_0, y_min = prop_y_min)
        t_deviation_1 = 1 - np_get_fraction(np.nansum(t_mc_permeths[t_demeth_ix], axis = 0), t_mc_total_1, y_min = prop_y_min)

        ### Do weighted average for the positions
        t_pop_mc_count = np.array(mc_count[:,ef_pop[1]], dtype = float).copy()
        t_pop_mc_total = np.array(mc_total[:,ef_pop[1]], dtype = float).copy()
        t_pop_mc_count[t_pop_mc_total <= mc_total_min] = np.nan
        t_pop_mc_total[t_pop_mc_total <= mc_total_min] = np.nan
        t_mc_read_total_0 = np.nansum( t_pop_mc_total[t_denovo_ix], axis = 0 )
        t_mc_read_total_1 = np.nansum( t_pop_mc_total[t_demeth_ix], axis = 0 )
        t_read_deviation_0 = np_get_fraction(np.nansum( t_pop_mc_count[t_denovo_ix], axis = 0 ), t_mc_read_total_0, y_min = prop_y_min)
        t_read_deviation_1 = 1 - np_get_fraction(np.nansum( t_pop_mc_count[t_demeth_ix], axis = 0 ), t_mc_read_total_1, y_min = prop_y_min)
        
        epi_out['deviations'] = epi_out['deviations'].append(
            pd.DataFrame( { 
                'subpop': ef_pop[0], 
                'deviation_0': t_read_deviation_0, 
                'deviation_1': t_read_deviation_1,
                'mc_total_0': t_mc_read_total_0,
                'mc_total_1': t_mc_read_total_1,
                'site_deviation_0': t_deviation_0, 
                'site_deviation_1': t_deviation_1,
                'site_total_0': t_mc_total_0,
                'site_total_1': t_mc_total_1
            },
            index = np.array(ef_pop[1]) )
        )
    return( epi_out )


def cache_pd_variable_to_file(var_key, cache_file, input_variable = None, mode = 'a'):
    if cache_file is not None:
        if op.exists(str(cache_file)):
            if var_key in h5.File(cache_file, 'r').keys():
                write_to_cache = False
                return( pd.read_hdf( cache_file, var_key )  )
    if input_variable is not None:
        input_variable.to_hdf( cache_file, key = var_key, mode = mode )
    return(input_variable)




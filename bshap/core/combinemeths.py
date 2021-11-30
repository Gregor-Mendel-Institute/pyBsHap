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

from . import meth5py
from pygenome import genome as g

def get_fraction(x, y, y_min = 0):
    if y <= y_min:
        return(np.nan)
    return(float(x)/y)

np_get_fraction = np.vectorize(get_fraction, excluded = "y_min")

class CombinedMethsTable(object):
    ## class for a combined meths tables

    def __init__(self, file_paths, file_ids = None, genome = "at_tair10"):
        self.meths_list = self.load_meths_files(file_paths)
        self.num_lines = len(self.meths_list)
        self.genome = g.GenomeClass(genome)
        if file_ids is not None:
            self.file_ids = file_ids
        else:
            self.file_ids = file_paths

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

    def derive_most_common_positions_echr(self, chrid, num_lines_with_info = 2):
        # returns common positions for a single chromosome
        #from tempfile import mkdtemp
        chrid_ind = self.genome.get_chr_ind(chrid)
        common_positions = np.arange(1, self.genome.golden_chrlen[chrid_ind] + 1, dtype="int")
        common_positions_scores = np.zeros( len(common_positions), dtype="int16" )
        t_echr_pos_ix = np.repeat(-1, len(common_positions) * self.num_lines).reshape((len(common_positions), self.num_lines))
        for m_ix in range(self.num_lines):
            m = self.meths_list[m_ix]
            e_chrinds = m.get_chrinds(chrid)
            e_chr_pos = m.positions[e_chrinds[1][0]:e_chrinds[1][1]]  ## these are the indices here
            common_positions_scores[e_chr_pos - 1] += 1
            t_echr_pos_ix[e_chr_pos - 1, m_ix] = np.arange(e_chrinds[1][0], e_chrinds[1][1])
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
            e_common_pos = self.derive_most_common_positions_echr(e_chr, num_lines_with_info)
            common_chrs = np.append(common_chrs, np.repeat( e_chr, e_common_pos[0].shape[0] ) )
            common_positions = np.append(common_positions, e_common_pos[0] )
            filter_pos_ixs = np.append(filter_pos_ixs, e_common_pos[1], axis = 0)
        log.info("done! total number of positions: %s" % str(len(common_positions)))
        self.write_methylated_positions(common_chrs, common_positions, filter_pos_ixs, output_file, min_mc_total = min_mc_total)
        return( h5.File(output_file, 'r' ) )

    def write_methylated_positions(self, common_chrs, common_positions, filter_pos_ixs, output_file, min_mc_total = 3, chunk_size=100000):
        num_positions = common_positions.shape[0]
        chunk_size = min( chunk_size, num_positions )
        log.info("writing the data into h5 file")
        outh5file = h5.File(output_file, 'w')
        data_mc_class = np.repeat( 'nan', num_positions )
        outh5file.create_dataset('chunk_size', data=chunk_size, shape=(1,),dtype='i8')
        outh5file.create_dataset('num_positions', data=num_positions, shape=(1,),dtype='i4')
        outh5file.create_dataset('file_ids', data = np.array(self.file_ids, dtype = "S"))
        outh5file.create_dataset('chr', data = np.array(common_chrs, dtype = "S"))
        outh5file.create_dataset('start', dtype = int, data = common_positions)
        outh5file.create_dataset('end', dtype = int, data = common_positions + 1)
        outh5file.create_dataset('filter_pos_ix', dtype = int, data = filter_pos_ixs, chunks=((chunk_size, self.num_lines)))
        mcs_count = outh5file.create_dataset('mc_count', fillvalue = 0, shape = (num_positions, self.num_lines), dtype = int, chunks=((chunk_size, self.num_lines)))
        mcs_total = outh5file.create_dataset('mc_total', fillvalue = 0, shape = (num_positions, self.num_lines), dtype = int, chunks=((chunk_size, self.num_lines)))
        mcs_permeths = outh5file.create_dataset('permeths', fillvalue = np.nan, shape = (num_positions, self.num_lines), dtype = float, chunks=((chunk_size, self.num_lines)))
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

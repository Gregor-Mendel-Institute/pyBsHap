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
from . import genome as g
genome = g.ArabidopsisGenome()


class CombinedMethsTable(object):
    ## class for a combined meths tables

    def __init__(self, file_paths):
        self.meths_list = self.load_meths_files(file_paths)
        self.num_lines = len(self.meths_list)

    def load_meths_files(self, file_paths):
        log.info("reading input files")
        meths_list = []
        for e_path in file_paths:
            meths_list.append(meth5py.load_hdf5_methylation_file(e_path))
        log.info("done! input %s files..!" % len(meths_list))
        return(meths_list)

    def derive_most_common_positions_echr(self, chrid, perc_nas):
        # returns common positions for a single chromosome
        #from tempfile import mkdtemp
        chrid_ind = genome.get_chr_ind(chrid)
        echr_pos_ix_file = op.join("chr_%s_inds" % str(chrid_ind + 1))
        if op.isfile(echr_pos_ix_file + '.npz'):
            echr_pos_ix = np.load( echr_pos_ix_file + '.npz' )
            return((echr_pos_ix['common_positions'], echr_pos_ix_file))
        common_positions = np.arange(1, genome.golden_chrlen[chrid_ind] + 1, dtype="int")
        common_positions_scores = np.zeros( len(common_positions), dtype="int16" )
        t_echr_pos_ix = np.repeat(-1, len(common_positions) * self.num_lines).reshape((len(common_positions), self.num_lines))
        for m_ix in range(self.num_lines):
            m = self.meths_list[m_ix]
            e_chrinds = m.get_chrinds(chrid)
            e_chr_pos = m.positions[e_chrinds[1][0]:e_chrinds[1][1]]  ## these are the indices here
            common_positions_scores[e_chr_pos - 1] += 1
            t_echr_pos_ix[e_chr_pos - 1, m_ix] = np.arange(len(e_chr_pos))
        ## Now calculate which positions are needed
        req_pos_ix = np.where( common_positions_scores >= perc_nas * self.num_lines )[0]
        common_positions = common_positions[req_pos_ix]
        np.savez( echr_pos_ix_file, common_positions = common_positions, pos_ix = t_echr_pos_ix[req_pos_ix, :] )
        return((common_positions, echr_pos_ix_file ))

    def derive_most_common_positions(self, perc_nas = 0.05):
        ## derive most common positions for all the chromosomes
        self.common_chrs_inds = []
        self.common_positions = np.zeros(0, dtype = "int" )
        self.file_pos_ix = []
        for e_chr in genome.chrs:
            log.info("reading in through chromosome %s" % e_chr)
            e_common_pos = self.derive_most_common_positions_echr(e_chr, perc_nas)
            self.common_chrs_inds.append( [len(self.common_positions), len(self.common_positions) + len(e_common_pos[0]) ] )
            self.common_positions = np.append( self.common_positions, e_common_pos[0] )
            self.file_pos_ix.append(e_common_pos[1])
        log.info("done! total number of positions: %s" % str(len(self.common_positions)))

    def get_methylated_values_pos_ix(self, filter_pos_inds):
        t_m_values = np.repeat(-1, filter_pos_inds.shape[0] * filter_pos_inds.shape[1]).astype("int8").reshape(filter_pos_inds.shape)
        for m_ix in range(self.num_lines):
            m = self.meths_list[m_ix]
            t_nan_ix = np.where( filter_pos_inds[:,m_ix] != -1 )[0]
            if len(t_nan_ix) > 0:
                t_m_values[t_nan_ix, m_ix] = m.__getattr__('methylated', filter_pos_inds[t_nan_ix,m_ix])
        return(t_m_values)

    def write_methylated_positions(self, output_file, read_depth=0, chunk_size=1000):
        ### Check if
        if 'common_positions' not in self.__dict__:
            log.error("please identify common_positions before writing the data into file")
        num_positions = len(self.common_positions)
        file_names = np.array([ op.basename(m.h5file.filename.encode('utf8')) for m in self.meths_list ], dtype="string")
        common_chrs = np.zeros( num_positions, dtype="S8" )
        for eind in range(len(genome.chrs)):
            common_chrs[ self.common_chrs_inds[eind][0]:self.common_chrs_inds[eind][1] ] = genome.chrs[eind]
        log.info("writing the data into h5 file")
        outh5file = h5.File(output_file, 'w')
        outh5file.create_dataset('chunk_size', data=chunk_size, shape=(1,),dtype='i8')
        outh5file.create_dataset('num_positions', data=num_positions, shape=(1,),dtype='i4')
        outh5file.create_dataset('files', data=file_names, shape=(self.num_lines,))
        outh5file.create_dataset('accessions', shape= ((self.num_lines,)), dtype='S20')
        outh5file.create_dataset('chr', compression="lzf", shape=(num_positions,), data = common_chrs)
        outh5file.create_dataset('start', compression="lzf", shape=(num_positions,), data = self.common_positions)
        outh5file.create_dataset('end', compression="lzf", shape=(num_positions,), data = self.common_positions + 1)
        meth_values = outh5file.create_dataset('value', compression = "lzf", shape = ((num_positions, self.num_lines)), chunks=((chunk_size, self.num_lines)),dtype='int')
        for ef, eg_ix in itertools.izip(self.file_pos_ix, self.common_chrs_inds):
            t_echr_pos_ix = np.load(ef + '.npz')['pos_ix']
            for ef_ix in range(0, eg_ix[1], chunk_size):
                t_m_values = self.get_methylated_values_pos_ix(t_echr_pos_ix[ef_ix:ef_ix+chunk_size,:])
                meth_values[eg_ix[0]+ef_ix:eg_ix[0]+ef_ix+t_m_values.shape[0],:] = t_m_values
            log.info("written data from file: %s" % ef)
        for m_ix in range(self.num_lines):
            outh5file['accessions'][m_ix] = file_names[m_ix].replace( "allc_", "" ).replace(".hdf5", "")
            ## generally id is like this allc_SRR3299777.hdf5
        outh5file.close()
        log.info("done")

    def derive_methylated_identical_pos_ix_echr(self, chrid, read_depth, pos_in_atleast, max_read_depth = 80):
        ## Caution: would probably need much memory if you give many files.
        # meth_list is a list of meths meth5py object read
        # returns common positions.
        chrid_ind = genome.get_chr_ind(chrid)
        num_positions = genome.golden_chrlen[chrid_ind]
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
        for e_chr in genome.chrs:
            # the function below is appending the list m.filter_pos_ix
            # So, make sure you have some place in the RAM.
            log.info("reading in through chromosome %s" % e_chr)
            self.derive_methylated_identical_pos_ix_echr(e_chr, read_depth=read_depth, pos_in_atleast=pos_in_atleast, max_read_depth=max_read_depth)
        log.info("done!")

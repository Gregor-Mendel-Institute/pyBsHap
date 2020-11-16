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
from . import genome as g
import csv
import itertools
import re

log = logging.getLogger(__name__)


def die(msg):
  sys.stderr.write('Error: ' + msg + '\n')
  sys.exit(1)

def get_chrs(allc_id, allc_path):
    #files_present = glob.glob(allc_path + "/allc_" + allc_id + "*tsv")
    if op.isfile(allc_path + "/allc_" + allc_id + "_Chr1.tsv"):
        return(['Chr1','Chr2','Chr3','Chr4','Chr5'])
    elif op.isfile(allc_path + "/allc_" + allc_id + "_1.tsv"):
        return(['1','2','3','4','5'])
    else:
        die("Either the allc id, %s or allc file path, %s is wrong" % (allc_id, allc_path))

def read_allc_pandas_table(allc_file):
    sniffer = csv.Sniffer()
    if op.splitext(allc_file)[1] == '.gz':  ### For the new allc files
        bsbed = pd.read_table(allc_file, header = None, compression='gzip')
    else:
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
        raise(NotImplementedError)
    return(bsbed)

def read_allc_files_chrwise(allc_id, allc_path):
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
    return((allcBed, chrpositions))

def generage_h5file_from_allc(allc_id, allc_path, outFile, ref_genome = "at_tair10", pval_thres=0.01):
    genome = g.ArabidopsisGenome(ref_genome)
    if allc_path == 'new':
        log.info("reading the allc file")
        allcBed = read_allc_pandas_table(allc_id)
        log.info("done!")
        if len(np.unique(allcBed.iloc[0:100000,6])) == 1: ### Check if the top one million is all either 1's or 0's
            from . import bsseq
            pval = bsseq.callMPs_allcbed(allcBed)
            allcBed.iloc[np.where(pval > pval_thres)[0],6] = 0   ### changing the methylated column to 0 for non methylated sites.
        allcBed_new = pd.DataFrame(columns = allcBed.columns)
        chrpositions = np.zeros(1, dtype="int")
        log.info("sorting the bed file!")
        sample_chrs = np.sort(np.unique(allcBed['chr']).astype("string"))
        sample_chrs_ix = [ genome.get_chr_ind(ec) for ec in sample_chrs ]
        if len(np.where(np.isfinite(np.array(sample_chrs_ix, dtype="float")))[0]) != len(genome.chrs):
            log.warn("The chromosome IDs do not match the tair IDs. Please check the file if they are suppposed to.")
        for ec in sample_chrs:
            allcBed_echr = allcBed.iloc[np.where(allcBed.iloc[:,0] == ec)[0], :]
            allcBed_new = allcBed_new.append(allcBed_echr.iloc[np.argsort(np.array(allcBed_echr.iloc[:,1], dtype=int)), :] , ignore_index=True)
            chrpositions = np.append(chrpositions, chrpositions[-1] + allcBed_echr.shape[0])
        log.info("done!")
        allcBed = allcBed_new
    else:
        (allcBed, chrpositions) = read_allc_files_chrwise(allc_id, allc_path)
    log.info("writing a hdf5 file")
    generate_H5File(allcBed,chrpositions, outFile)

def generate_H5File(allcBed, chrpositions, outFile):
    chunk_size = min(100000, allcBed.shape[0])
    h5file = h5.File(outFile, 'w')
    h5file.create_dataset('chunk_size', data=chunk_size, shape=(1,),dtype='i8')
    h5file.create_dataset('chrpositions', data=chrpositions, shape=(len(chrpositions),),dtype='int')
    allc_columns = ['chr', 'pos', 'strand', 'mc_class', 'mc_count', 'total', 'methylated']
    allc_dtypes = ['S16', 'int', 'S4', 'S8', 'int', 'int', 'int8']
    ## Going through all the columns
    for cols, coltype in zip(allc_columns, allc_dtypes):
        h5file.create_dataset(cols, compression="gzip", data=np.array(allcBed[cols],dtype=coltype), shape=(allcBed.shape[0],), chunks = ((chunk_size,)))
    if allcBed.shape[1] > 7:
        for i in range(7,allcBed.shape[1]):
            try:
                h5file.create_dataset(allcBed.columns[i], compression="gzip", data=np.array(allcBed[allcBed.columns[i]]), shape=(allcBed.shape[0],),chunks = ((chunk_size,)))
            except TypeError:
                h5file.create_dataset(allcBed.columns[i], compression="gzip", data=np.array(allcBed[allcBed.columns[i]],dtype="S8"), shape=(allcBed.shape[0],),chunks = ((chunk_size,)))
    h5file.close()



def load_hdf5_methylation_file(hdf5_file, ref_genome="at_tair10", bin_bed=''):
    return HDF5MethTable(hdf5_file, ref_genome, bin_bed)

# Try to make it as a class, learned from PyGWAS
class HDF5MethTable(object):

    def __init__(self,hdf5_file, ref_genome = "at_tair10", bin_bed=''):
        self.h5file = h5.File(hdf5_file,'r')
        self.filter_pos_ix = self.get_filter_inds(bin_bed)
        self.chrpositions = np.array(self.h5file['chrpositions'])
        self.chrs = self.__getattr__('chr', self.chrpositions[0:-1])
        self.genome =  g.ArabidopsisGenome(ref_genome)

    def close_h5file(self):
        if self.h5file is not None:
            self.h5file.close()

    def get_chrinds(self, chrid):
        chrid_mod = str(chrid).replace("Chr", "").replace("chr", "")
        real_chrs = np.array( [ ec.replace("Chr", "").replace("chr", "") for ec in self.chrs ] )
        req_chr_ind = np.where(real_chrs == chrid_mod)[0]
        if len( req_chr_ind ) == 0:
            return( (None, None) )
        req_chr_ind = req_chr_ind[0]
        chr_inds = [self.chrpositions[req_chr_ind], self.chrpositions[req_chr_ind + 1]]
        return((req_chr_ind, chr_inds))

    def get_filter_inds(self, bin_bed):
        # bin_bed = ['Chr1', 0, 100] or "Chr1,1,100"
        if bin_bed == '':
            return(None)
        if isinstance(bin_bed, basestring):
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
                return(np.array(self.h5file[str(name)]))
            return(self.h5file[str(name)])
        elif type(filter_pos_ix) is np.ndarray:
            if len(filter_pos_ix) == 0:
                ret_attr = np.array(())
            else:
                rel_pos_ix = filter_pos_ix - filter_pos_ix[0]
                ret_attr = np.array(self.h5file[str(name)][filter_pos_ix[0]:filter_pos_ix[-1]+1][rel_pos_ix])
        else:
            ret_attr = np.array(self.h5file[str(name)][filter_pos_ix])
        if name in ['chr', 'mc_class', 'strand']:
            ret_attr = ret_attr.astype('U13')
        elif name in ['pos', 'total', 'mc_count', 'methylated']:
            ret_attr = ret_attr.astype(int)
        return(ret_attr)

    def get_bed_df(self, filter_pos_ix, full_bed=False, read_threshold=0):
        req_chrs = self.__getattr__('chr', filter_pos_ix, return_np=True)
        req_pos = self.__getattr__('pos', filter_pos_ix, return_np=True)
        if full_bed:
            conname = self.__getattr__('mc_class', filter_pos_ix, return_np=True )
            strand = self.__getattr__('strand', filter_pos_ix, return_np=True)
            permeth = self.get_permeths(filter_pos_ix, read_threshold=read_threshold)
            return(pd.DataFrame(np.column_stack((req_chrs, req_pos, req_pos + 1, conname, permeth, strand)), columns=['chr', 'start', 'end', 'mc_class', 'permeth', 'strand']))
        return(pd.DataFrame(np.column_stack((req_chrs, req_pos, req_pos + 1)), columns=['chr', 'start', 'end']))

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
        methylated_cs = np.array(self.__getattr__('methylated', filter_pos_ix), dtype="float")
        mc_count = self.__getattr__('mc_count', filter_pos_ix, return_np=True)
        mc_total = self.__getattr__('mc_total', filter_pos_ix, return_np=True)
        calculated_permeths = np.divide(np.multiply(methylated_cs, mc_count), mc_total)
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
            bin_bed = [int(t), int(t) + window_size - 1]
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

    def MethylationSummaryStats(self, filter_pos_ix, category = 1, req_context = None):
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
            methylated_cs = np.where(self.__getattr__('methylated', filter_pos_ix, return_np=True) == 1)[0]
            return(np.divide(float(np.sum(mc_count[methylated_cs])), np.sum(mc_total)))
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
                if count % 1000 == 0:
                    log.info("progress: analysed %s windows" % count)
        outmeths_cg_avg.close()
        outmeths_chg_avg.close()
        outmeths_chh_avg.close()
        log.info("done!")

    def generate_meth_average_required_bed(self, required_bed, out_file, category=1):
        if op.isfile(required_bed):
            req_regions = pd.read_csv(required_bed, sep = "\t", header=None)
        elif type(required_bed) is str:
            ### required_bed = "Chr1,1,1000"
            req_regions = pd.DataFrame(required_bed.split(",")).T
        elif type(required_bed) is pd.DataFrame:
            req_regions = required_bed
        assert len(req_regions.columns) >= 3, "bed file should have minimum of 3 columns. eg., Chr, Start, End"
        outmeths_cg_avg = open(out_file + ".CG.bg", 'w')
        outmeths_chg_avg = open(out_file + ".CHG.bg", 'w')
        outmeths_chh_avg = open(out_file + ".CHH.bg", 'w')
        for er in req_regions.iterrows():
            t_filter_pos_ix = self.get_filter_inds( [er[1][0], int(er[1][1]), int(er[1][2]) ] )
            count = 0
            t_count = self._TotalCounts_All_Contexts( t_filter_pos_ix )
            req_meth_avg = self._AveMethylation_All_Contexts(t_filter_pos_ix, category)
            outmeths_cg_avg.write("%s\t%s\t%s\t%s\t%s\n" % (er[1][0], int(er[1][1]), int(er[1][2]), req_meth_avg[0], len(t_count[0])))
            outmeths_chg_avg.write("%s\t%s\t%s\t%s\t%s\n" % (er[1][0], int(er[1][1]), int(er[1][2]), req_meth_avg[1], len(t_count[1])))
            outmeths_chh_avg.write("%s\t%s\t%s\t%s\t%s\n" % (er[1][0], int(er[1][1]), int(er[1][2]), req_meth_avg[2], len(t_count[2])))
            count = count + 1
            if count % 100 == 0:
                log.info("progress: analysed %s regions" % count)
        outmeths_cg_avg.close()
        outmeths_chg_avg.close()
        outmeths_chh_avg.close()
        log.info("done!")

def load_input_meths(inFile, allc_path, ref_genome):
    assert op.isfile(inFile), "input file is not present"
    _,inType = op.splitext(inFile)
    if inType == '.hdf5':
        log.info("reading hdf5 file!")
        meths = load_hdf5_methylation_file(inFile, ref_genome)
        log.info("done!")
        return(meths)
    outhdf5 = op.splitext(op.splitext(inFile)[0])[0] + ".hdf5"
    if len(re.compile(".tsv.gz$").findall(inFile)) > 0:
        log.info("generating hdf5 file from allc files, %s!" % outhdf5)
        generage_h5file_from_allc(inFile, "new", outhdf5)
    log.info("loading the hdf5 file %s!" % outhdf5)
    meths = load_hdf5_methylation_file(outhdf5, ref_genome)
    log.info("done!")
    return(meths)

def potatoskin_methylation_averages(args):
    # bin_bed = Chr1,1,100
    if args['allc_path'] == "bed":
        args['outFile'] =  args['outFile'] + '.bedGraph'
        run_bedtools.get_genomewide_methylation_average(args)
        log.info("finished!")
        return(0)
    meths = load_input_meths(args['inFile'], args['allc_path'], args['ref_genome'])
    if args['required_region'] == '0,0,0':
        meths.generate_meth_average_in_windows(args['outFile'], args['window_size'], category=args['category'])
        return(0)
    meths.generate_meth_average_required_bed(args['required_region'], args['outFile'], args['category'])
    log.info('finished!')
    return(0)

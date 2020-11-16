# Loading genome data
import logging
import numpy as np
import pandas as pd
import string
import sys, os.path

log = logging.getLogger(__name__)

def die(msg):
  sys.stderr.write('Error: ' + msg + '\n')
  sys.exit(1)

def getInd_bin_bed(bin_bed, tair10):
    ## bin_bed = ["Chr1", 1, 1000]
    bin_s = [int(bin_bed[0].replace("Chr", "")) - 1, int(bin_bed[1]), int(bin_bed[2])]
    return(tair10.chr_inds[bin_s[0]] + int(( bin_s[1] + bin_s[2] )/2) )

class ArabidopsisGenome(object):
    ## coordinates for ArabidopsisGenome using TAIR 10

    def __init__(self, ref_genome = "at_tair10"):
        if ref_genome == "at_tair10":
            self.chrs = ['Chr1','Chr2','Chr3','Chr4','Chr5']
            self.def_color = ["#1f78b4", "#33a02c", "#1f78b4", "#33a02c", "#1f78b4"]
            self.real_chrlen = [34964571, 22037565, 25499034, 20862711, 31270811]
            self.golden_chrlen = [30427671, 19698289, 23459830, 18585056, 26975502]
            self.centro_start = [14364752, 3602775, 12674550, 2919690, 11668616]
            self.centro_end   = [15750321, 3735247, 13674767, 4011692, 12082583]
            self.cetro_mid = np.add(self.centro_start, self.centro_end)/2
        elif os.path.exists(ref_genome):
            ## Provide a fasta file to check for genome lengths etc
            from pyfaidx import Faidx
            genome = Faidx(ref_genome).index
            self.chrs = np.sort(np.array([ef for ef in genome.keys()])).tolist()
            self.real_chrlen = [ genome[ef].rlen for ef in self.chrs]
            self.golden_chrlen = self.real_chrlen
        self.chr_inds = np.append(0, np.cumsum(self.golden_chrlen))

    def load_genome_fasta(self, fasta_file):
        from pyfaidx import Fasta
        self.fasta = Fasta(fasta_file)

    def load_bed_ids_str(self, **kwargs):
        for req_name in kwargs:
            req_bed_df = pd.read_csv( kwargs[req_name], header=None, sep = "\t")
            setattr(self, req_name, req_bed_df)
            setattr(self, req_name + "_str", np.array(req_bed_df.iloc[:,0] + ',' + req_bed_df.iloc[:,1].map(str) + ',' + req_bed_df.iloc[:,2].map(str), dtype="str") )

    def determine_bed_from_araportids(self, name, araportids, return_bed=False):
        assert type(araportids) is pd.core.series.Series, "please provide a pandas series object"
        assert hasattr(self, name), "please load required bed file using 'load_bed_ids_str' function. ex., ARAPORT11/Araport11_GFF3_genes_201606.bed"
        bed_str = np.zeros(0, dtype="float")
        bed_ix = np.zeros(0, dtype="int")
        for ei in araportids:
            t_ind = np.where( self.__getattribute__( name ).iloc[:,3] == ei )[0]
            if len(t_ind) == 0:
                bed_str = np.append(bed_str, '')
            else:
                bed_str = np.append( bed_str, self.__getattribute__( name + "_str" )[t_ind[0]] )
                bed_ix = np.append(bed_ix, t_ind[0])
        if return_bed:
            bed_str = self.__getattribute__( name ).iloc[bed_ix,]
            bed_str = bed_str.rename(columns={0:"chrom", 1:"start", 2:"end", 3:'id', 4:'score',5:'strand'})
        return( bed_str )

    def get_chr_ind(self, echr):
        real_chrs = np.array( [ ec.replace("Chr", "").replace("chr", "") for ec in self.chrs ] )
        if type(echr) is str or type(echr) is np.string_:
            echr_num = str(echr).replace("Chr", "").replace("chr", "")
            if len(np.where(real_chrs == echr_num)[0]) == 1:
                return(np.where(real_chrs == echr_num)[0][0])
            else:
                return(None)
        echr_num = np.unique( np.array( echr ) )
        ret_echr_ix = np.zeros( len(echr), dtype="int8" )
        for ec in echr_num:
            t_ix = np.where(real_chrs ==  str(ec).replace("Chr", "").replace("chr", "") )[0]
            ret_echr_ix[ np.where(np.array( echr ) == ec)[0] ] = t_ix[0]
        return(ret_echr_ix)

    def get_genomewide_inds(self, df_str):
        ### This is the function to give the indices of the genome when you give a bed file.
        if type(df_str) is not pd.core.series.Series and type(df_str) is not pd.core.frame.DataFrame:
            die("please input pandas dataframe or series object")
        elif type(df_str) is pd.core.series.Series:
            df_str_np = np.array(df_str, dtype="string")
            df_str_unique = np.unique(df_str_np, return_inverse=True)
            df_str_inds = np.array(pd.Series(df_str_unique[0]).str.split(",").apply(getInd_bin_bed, args= (self,) ))
            return( df_str_inds[df_str_unique[1]] )
        elif type(df_str) is pd.core.frame.DataFrame:
            ## here first column is chr and second is position
            if df_str.shape[1] == 3:
                f_df_str = pd.DataFrame(df_str.iloc[:,0])
                f_df_str['pos'] = pd.Series( ((df_str.iloc[:,1] + df_str.iloc[:,2]) / 2).apply(int) )
            else:
                f_df_str = df_str
            chrom_ix = np.array(self.chr_inds[self.get_chr_ind( f_df_str.iloc[:,0] )], dtype=int)
            return( chrom_ix + np.array(f_df_str.iloc[:,1], dtype=int) )

    def get_genomic_position_from_ind(self, ind):
        ## This function is just opposite to the one before.
        # Given an index, we should get the chromosomes and position
        chr_idx = np.searchsorted( self.chr_inds, ind ) - 1
        ind_pos = ind - self.chr_inds[chr_idx]
        ind_chr = np.array(self.chrs)[chr_idx]
        return( pd.DataFrame( np.column_stack((ind_chr, ind_pos)), columns = ["chr", "pos"]  ) )

    def iter_windows_echr(self, echr, window_size, overlap=0):
        chr_ind = self.get_chr_ind(echr)
        if overlap >= window_size:
            raise(NotImplementedError)
        if overlap > 0:
            for x in range(1, self.golden_chrlen[chr_ind], overlap):
                yield([x, x + window_size - 1])
        else:
            for x in range(1, self.golden_chrlen[chr_ind], window_size):
                yield([x, x + window_size - 1])

    def iter_windows(self, window_size, overlap=0):
        for echr, echrlen in zip(self.chrs, self.golden_chrlen):
            echr_windows = self.iter_windows_echr(echr, window_size, overlap)
            for ewin in echr_windows:
                yield([echr, ewin[0], ewin[1]])

    def estimated_cM_distance(self, snp_position):
        ## snp_position = "Chr1,150000" or "Chr1,1,300000"
        # Data based on
        #Salome, P. A., Bomblies, K., Fitz, J., Laitinen, R. A., Warthmann, N., Yant, L., & Weigel, D. (2011)
        #The recombination landscape in Arabidopsis thaliana F2 populations. Heredity, 108(4), 447-55.
        assert isinstance(snp_position, basestring)
        assert len(snp_position.split(",")) >= 2
        if len(snp_position.split(",")) == 2:
            snp_position = [snp_position.split(",")[0], int(snp_position.split(",")[1])]
        elif len(snp_position.split(",")) == 3:
            snp_position = [snp_position.split(",")[0], (int(snp_position.split(",")[1]) + int(snp_position.split(",")[2])) / 2 ]
        mean_recomb_rates = [3.4, 3.6, 3.5, 3.8, 3.6]  ## cM/Mb
        chr_ix = self.get_chr_ind( snp_position[0] )
        return( mean_recomb_rates[chr_ix] * snp_position[1] / 1000000 )

    def get_mc_context(self, cid, pos):
        dnastring_pos = self.fasta[cid][pos:pos+3].seq.encode('ascii').upper()
        ## make sure you can have to identify strand here
        dnastring_neg = self.fasta[cid][pos-2:pos+1].seq.encode('ascii').upper()  ## Changed the position, different from forward
        dnastring_neg = get_reverse_complement(dnastring_neg)
        if dnastring_pos[0].upper() == 'C' and dnastring_neg[0].upper() != 'C':
            strand = '0'
            dnastring = dnastring_pos
        elif dnastring_pos[0].upper() != 'C' and dnastring_neg[0].upper() == 'C':
            strand = '1'
            dnastring = dnastring_neg
        else:
            return((dnastring_pos, dnastring_neg, None))
        if dnastring[1].upper() == 'G':
            dna_context = ["CG",0]
        elif dnastring[2].upper() == 'G':
            dna_context = ["CHG",1]
        elif dnastring:
            dna_context = ["CHH",2]
        return((dnastring, dna_context, strand))

    def get_inds_overlap_region(self, region_bed_df, name="genes", request_ind = None, g = None):
        import pybedtools as pybed
        assert hasattr(self, name), "please load required bed file using 'load_bed_ids_str' function. ex., ARAPORT11/Araport11_GFF3_genes_201606.bed"
        assert type(region_bed_df) is pd.core.frame.DataFrame, "please provide a pandas series object"
        if request_ind is None:
            gene_bed = pybed.BedTool.from_dataframe( self.__getattribute__(name).iloc[:,[0,1,2]] )
        else:
            gene_bed = pybed.BedTool.from_dataframe( self.__getattribute__(name).iloc[:,[0,1,2]].iloc[[request_ind]]  )
        region_bed = pybed.BedTool.from_dataframe(region_bed_df.iloc[:,[0,1,2]])
        region_bed_str = np.array(region_bed_df.iloc[:,0].map(str) + "," + region_bed_df.iloc[:,1].map(str) + "," +  region_bed_df.iloc[:,2].map(str), dtype="str")
        ## Just taking first three columns for bedtools
        inter_region_bed = region_bed.intersect(gene_bed, wa=True)
        inter_gene_bed = gene_bed.intersect(region_bed, wa=True)
        if inter_region_bed.count() == 0:   ## Return if there are no matching lines.
            return(None)
        inter_region_bed = inter_region_bed.to_dataframe() ## wa is to return the entire bed.
        inter_region_bed_str = np.array(inter_region_bed.iloc[:,0].map(str) + "," + inter_region_bed.iloc[:,1].map(str) + "," +  inter_region_bed.iloc[:,2].map(str), dtype="str")
        inter_gene_bed = inter_gene_bed.to_dataframe()
        inter_gene_bed_str = np.array(inter_gene_bed.iloc[:,0].map(str) + "," + inter_gene_bed.iloc[:,1].map(str) + "," +  inter_gene_bed.iloc[:,2].map(str), dtype="str")
        out_dict = { "region_ix": np.where( np.in1d(region_bed_str, inter_region_bed_str ) )[0] }
        out_dict['ref_ix'] = np.where( np.in1d(self.__getattribute__(name + '_str'), inter_gene_bed_str ) )[0]
        return(out_dict)

def get_reverse_complement(seq):
    old_chars = "ACGT"
    replace_chars = "TGCA"
    tab = string.maketrans(old_chars,replace_chars)
    return(seq.translate(tab)[::-1])

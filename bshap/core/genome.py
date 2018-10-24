# Loading genome data
import logging
import numpy as np
import pandas as pd
import string


log = logging.getLogger(__name__)

def getInd_bin_bed(bin_bed, tair10):
    ## bin_bed = ["Chr1", 1, 1000]
    bin_s = [int(bin_bed[0].replace("Chr", "")) - 1, int(bin_bed[1]), int(bin_bed[2])]
    return(tair10.chr_inds[bin_s[0]] + int(( bin_s[1] + bin_s[2] )/2) )

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

    def load_genome_fasta(self, fasta_file):
        from pyfaidx import Fasta
        self.fasta = Fasta(fasta_file)

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
        except IndexError as err_idx:
            return(None)

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
                df_str = pd.DataFrame(df_str.iloc[:,0]).join(pd.DataFrame( ((df_str.iloc[:,1] + df_str.iloc[:,2]) / 2).apply(int) ))
            chrom = np.char.replace(np.core.defchararray.lower(np.array(df_str.iloc[:,0], dtype="string")), "chr", "")
            return(self.chr_inds[np.array(chrom, dtype=int) - 1] + np.array(df_str.iloc[:,1]) )

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

def get_reverse_complement(seq):
    old_chars = "ACGT"
    replace_chars = "TGCA"
    tab = string.maketrans(old_chars,replace_chars)
    return(seq.translate(tab)[::-1])

## Main module for pyBsHap

import pysam
from pyfaidx import Fasta
import logging
import numpy as np
import re
import h5py
import string
import os.path, sys
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
import json
from . import meth5py

log = logging.getLogger(__name__)

def die(msg):
  sys.stderr.write('Error: ' + msg + '\n')
  sys.exit(1)

def getChrs(inBam):
    chrs = np.array([x['SN'] for x in inBam.header['SQ']])
    chrs1d = np.char.replace(np.core.defchararray.lower(np.array(chrs, dtype="string")), "chr", "")
    chrslen = np.array([x['LN'] for x in inBam.header['SQ']])
    reqchrs = chrs[np.where(np.char.isdigit(chrs1d))[0]].tolist()
    reqchrslen = chrslen[np.where(np.char.isdigit(chrs1d))[0]].tolist()
    binLen = np.zeros(0, dtype = int)
    for tread in inBam.fetch(reqchrs[0], 0, reqchrslen[0]):
        binLen = np.append(binLen, tread.infer_query_length())
        if len(binLen) >= 100:
            break
    #return (reqchrs, reqchrslen, int(np.nanmean(binLen)))
    return(chrs, chrslen, int(np.nanmean(binLen)))

def get_reverse_complement(seq):
    old_chars = "ACGT"
    replace_chars = "TGCA"
    tab = string.maketrans(old_chars,replace_chars)
    return(seq.translate(tab)[::-1])

def findIntersectbed(bin_bed, map_bed):
    # bin_bed = [290605, 290688]
    # map_bed = [290626,290711]
    bin_bed_arr = np.array(list(range(bin_bed[0], bin_bed[1])))
    map_bed_arr = np.array(list(range(map_bed[0], map_bed[1])))
    intersect_bed_arr = np.sort(np.intersect1d(bin_bed_arr, map_bed_arr))
    intersect_bed = [intersect_bed_arr[0], intersect_bed_arr[-1]+1]
    return(intersect_bed)

def parseContext(tair10, cid, pos, strand):
    try:
        if strand == '0':
            dnastring = tair10[cid][pos:pos+3].seq.encode('ascii').upper()
        elif strand == '1': ## make sure you can have to identify strand here
            dnastring = tair10[cid][pos-2:pos+1].seq.encode('ascii').upper()  ## Changed the position, different from forward
            dnastring = get_reverse_complement(dnastring)
    except:
        return("CNN", "CN")
    if dnastring[1].upper() == 'G':
        dna_context = ["CG",0]
    elif dnastring[2].upper() == 'G':
        dna_context = ["CHG",1]
    elif dnastring:
        dna_context = ["CHH",2]
    return(dnastring, dna_context)

def decodeFlag(flag):
    bintag = str(int(bin(flag)[2:]))[::-1]
    notag = (12 - len(bintag)) * '0'
    return(bintag + notag)

def getStrandPE(binread):
    df = decodeFlag(binread.flag)
    if binread.is_read2:
        ## if it is read2, then we have to reverse strand flag? need to double check this one
        if df[4] == '1':
            return('0')
        else:
            return('1')
    else:
        return(df[4])

def getStrandContext(strand):
    if strand == '0': ## Forward read, count the number of C's
        re_string = ['C','T','.']
    elif strand == '1': ## reverse
        re_string = ['G', 'A',',']
    return(re_string)

def getHighlightedSeqs(refseq, rseq, strand):
    re_string = getStrandContext(strand)
    dot_rseq = rseq
    for i,c in enumerate(refseq):
        if c.upper() != re_string[0]:
            dot_rseq = dot_rseq[:i] + re_string[2] + dot_rseq[i+1:]
    return(dot_rseq)

def getMatchedSeq(binread, tair10, bin_bed):
    intersect_bed = findIntersectbed(bin_bed[1:3], [binread.reference_start, binread.reference_end])
    bin_positions = np.where(np.in1d(binread.get_reference_positions(), list(range(intersect_bed[0], intersect_bed[1]))))[0]
    refseq = ''
    rseq = ''
    for ind in bin_positions:
        (i, rind) = binread.get_aligned_pairs()[ind]
        if i is not None and rind is not None:
            rseq = rseq + binread.seq[i]
            refseq = refseq + tair10[bin_bed[0]][rind].seq.encode('ascii')
        elif i is None and rind is not None:
            rseq = rseq + '-'
            refseq = refseq + tair10[bin_bed[0]][rind].seq.encode('ascii')
        elif rind is None and i is not None:
            rseq = rseq + binread.seq[i]
            refseq = refseq + '-'
    return(refseq, rseq)

def getSeqRecord(binread, tair10, bin_bed, intersect_bed, defreturn = 0):
    ## bin_bed = ['Chr1', start, end, binLen, pointPos]
    ## intersect_bed = [s, e]
    char_add = 'N'
    strand = getStrandPE(binread)
    refseq, rseq = getMatchedSeq(binread, tair10, bin_bed)
    dot_rseq = getHighlightedSeqs(refseq, rseq, strand)
    fin_rseq = char_add * (intersect_bed[0] - bin_bed[1]) + dot_rseq + char_add * (bin_bed[2] - intersect_bed[1])
    if defreturn == 0:
        return(SeqRecord(Seq(fin_rseq, generic_dna), id = binread.query_name.split(' ')[0]))
    else:
        return(fin_rseq)

def getMethRead(tair10, binread):
    strand = getStrandPE(binread)
    error_rate = 0.01
    context = getStrandContext(strand)[0]
    mCs, tCs, read_length = [0,0,0], [0,0,0], 0
    for i,rind in binread.get_aligned_pairs():
        read_length = read_length + 1
        if i is None or rind is None:  ## Skip the position if there is an indel
            continue
        c = tair10[binread.reference_name][rind].seq.encode('ascii')
        if c.upper() == context:
            (dnastring, methcontext) = parseContext(tair10, binread.reference_name, rind, strand)
            tCs[methcontext[1]] += 1
            if binread.seq[i].upper() == context:
                mCs[methcontext[1]] += 1
    nc_thres = 3
    rmeth = [float(mCs[i])/tCs[i] if tCs[i] >= nc_thres else -1 for i in range(3)]
    tmeth = float(sum(mCs))/sum(tCs) if sum(tCs) >= nc_thres else -1
    #permeths = [tmeth, rmeth[0], rmeth[1], rmeth[2]]
    #permeths = [tmeth, tCs[0], tCs[1], tCs[2]]
    permeths = [tmeth, tCs[0], rmeth[0], tCs[1], rmeth[1], tCs[2], rmeth[2]]
    return(permeths)
    # We need to differentiate between the reads which are small and has no methylation

def filterRead(binread):
    rflag = decodeFlag(binread.flag)
    filter_flag = False
    if len(binread.get_aligned_pairs()) == 0 or rflag[10] != '0' or rflag[8] != '0':
        ## removing the duplicates
        ## Filter the read if alignment is not good!
        filter_flag = True
    return(filter_flag)

def getMethWind(inBam, tair10, required_bed, meths = '', strand = '0'):
    # required_bed = ['Chr1', start, end, binLen]
    binLen = required_bed[3]
    read_length_thres = binLen/2
    bin_start = required_bed[1] - (required_bed[1] % binLen)
    estimated_bins = list(range(bin_start, required_bed[2], binLen))
    #dt = h5py.special_dtype(vlen=np.dtype('float16'))
    #reqmeths = meths.create_dataset(required_bed[0],compression="gzip", shape = (len(estimated_bins),4,), dtype=dt)
    #meths[required_bed[0]].attrs['positions'] = estimated_bins
    if meths != '':
        meths.create_dataset("bins_" + required_bed[0], data = estimated_bins)
    progress_bins = 0
    #window_alignment = []  ## multiple segment alignment
    binmeth_whole = []
    for bins in estimated_bins:        ## sliding windows with binLen
        binmeth = []
        progress_bins += 1
        dot_refseq = tair10[required_bed[0]][bins:bins + binLen].seq.encode('ascii')
        dot_refseq = re.sub('A|T', '.', dot_refseq.upper())
        bins_alignment = [SeqRecord(Seq(dot_refseq, generic_dna), id = required_bed[0] + ':' + str(bins) + '-' + str(bins + binLen))]
        for binread in inBam.fetch(required_bed[0], bins, bins + binLen):
            r_strand = getStrandPE(binread)
            bin_bed = [required_bed[0], bins, bins + binLen]
            if filterRead(binread):
                continue
            intersect_bed = findIntersectbed(bin_bed[1:3], [binread.reference_start, binread.reference_end])
            intersect_len = len(list(range(intersect_bed[0], intersect_bed[1])))
            rseq_record = getSeqRecord(binread, tair10, bin_bed, intersect_bed) ## No need to check for the overlap to print the alignment
            if intersect_len > read_length_thres:  ## To print the read only in one window.
                permeths = getMethRead(tair10, binread) ## taking the first and last
                binmeth.append([permeths[i] for i in [0,2,4,6]])
                binmeth_whole.append(permeths)
                if r_strand == strand:         ### Get only forward reads in the alignment
                    bins_alignment.append(rseq_record)
        if progress_bins % 1000 == 0:
            log.info("ProgressMeter - %s windows in analysed, %s total" % (progress_bins, len(estimated_bins)))
        if np.array(binmeth).shape[0] > 0 and meths != '':
            meths.create_dataset("b_" + required_bed[0] + "_" + str(bins),compression="gzip", data = np.array(binmeth).T)
        #reqmeths[progress_bins-1] = np.array(binmeth).T
        #window_alignment.append(MultipleSeqAlignment(bins_alignment))
    return(binmeth_whole)

def getMethGenome(bamFile, fastaFile, outFile, interesting_region='0,0,0'):
    ## input a bam file, bin length
    ## make sure the binLen is less than the read length
    log.info("loading the input files!")
    inBam = pysam.AlignmentFile(bamFile, "rb")
    (chrs, chrslen, binLen) = getChrs(inBam)
    tair10 = Fasta(fastaFile)
    log.info("finished!")
    meths = h5py.File('meths.' + outFile + '.hdf5', 'w')
    meths.create_dataset('input_bam', data = os.path.basename(bamFile))
    meths.create_dataset('chrs', data = chrs)
    meths.create_dataset('binlen', data = binLen)
    meths.create_dataset('chrslen', data = chrslen)
    if interesting_region == '0,0,0':
        for cid, clen in zip(chrs, chrslen):     ## chromosome wise
            log.info("analysing chromosome: %s" % cid)
            required_bed = [cid, 1, clen, binLen]   ### 0 is to not print reads
            getMethWind(inBam, tair10, required_bed, meths)
            log.info("finished!")
            meths.close()
        return(0)
    else:
        required_region = interesting_region.split(',')
        required_bed = [required_region[0], int(required_region[1]), int(required_region[2]), binLen, 1]
        log.info("analysing region %s:%s-%s !" % (required_bed[0], required_bed[1], required_bed[2]))
        binmeth_whole = getMethWind(inBam, tair10, required_bed, meths)
        if len(required_region) == 3:
            #AlignIO.write(window_alignment, 'meths.' + outFile + '.aln', "clustal")
            (type_counts,type_freqs) = clusteringReads(binmeth_whole)
            types_summary = {'counts': type_counts, 'region': interesting_region, 'input_bam': os.path.basename(bamFile), 'freqs': type_freqs, 'binmeth_whole': binmeth_whole}
            with open('meths.' + outFile + '.summary.json', "w") as out_stats:
                out_stats.write(json.dumps(types_summary))
        meths.close()
        log.info("finished!")
        return 0

def getMethsRegions(bamFile, fastaFile, outFile, regionsFile):
    ## Give a bed file for this
    ## it should be tab delimited
    log.info("loading the input files!")
    inBam = pysam.AlignmentFile(bamFile, "rb")
    (chrs, chrslen, binLen) = getChrs(inBam)
    tair10 = Fasta(fastaFile)
    log.info("done!")
    outTxt = open('meths.' + outFile + '.summary.txt', 'w')
    with open(regionsFile) as rFile:
        for rline in rFile:
            rline_split = rline.rstrip().split('\t')
            required_bed = [rline_split[0], int(rline_split[1]), int(rline_split[2]), binLen]
            log.info("analysing region %s:%s-%s !" % (required_bed[0], required_bed[1], required_bed[2]))
            binmeth_whole = getMethWind(inBam, tair10, required_bed, '')
            (type_counts,type_freqs) = clusteringReads(binmeth_whole)
            outTxt.write("%s" % rline_split[0] + ',' + rline_split[1]  + ',' + rline_split[2])
            for i in range(len(type_counts)):
                outTxt.write(",%s" % type_counts[i])
            outTxt.write("\n")
    outTxt.close()
    log.info("finished!")

def getCmatrix(bins_alignment, strand):
    # re_string = ['C','T','.']
    re_string = getStrandContext(strand)
    num_seq = len(bins_alignment)
    refseq = bins_alignment[0]
    num_cs = str(refseq.seq).count(re_string[0])
    wind_hap = np.zeros((num_seq-1, num_cs), dtype="int8")
    wind_hap[wind_hap == 0] = -1
    ind = 0
    for ec in re.finditer(re_string[0], str(refseq.seq)):
        for seqind in range(1, num_seq):
            seqNu = str(bins_alignment[seqind].seq)[ec.start()]
            if seqNu.upper() == re_string[0]:
                wind_hap[seqind-1, ind] = 1
            elif seqNu.upper() == re_string[1]:
                wind_hap[seqind-1, ind] = 0
        ind = ind + 1
    return wind_hap

def haplotypeBlocks(bins_alignment, strand = '0'):
    ## Here strand is to check C's or G's
    num_seq = len(bins_alignment)
    if num_seq <= 1:
        return 0
    c_matrix = getCmatrix(bins_alignment, strand)
    return c_matrix

def countTypeFreqs(type_cls):
    numCLs = 8
    type_cls_freq = np.array(type_cls[1],dtype=float)/np.nansum(type_cls[1])
    type_counts = []
    type_freqs = []
    for i in range(numCLs):
        if len(np.where(type_cls[0] == i)[0]) > 0:
            type_counts.append(type_cls[1][np.where(type_cls[0] == i)[0][0]])
            type_freqs.append(type_cls_freq[np.where(type_cls[0] == i)[0][0]])
        else:
            type_counts.append(0)
            type_freqs.append(0)
    return type_counts,type_freqs

def clusteringReads(binmeth_whole, n_clusters=8):
    from sklearn.cluster import KMeans
    init_cls = np.array(((0,0,0),(1,0,0),(0,1,0), (0,0,1), (1,1,0), (1,0,1), (0,1,1), (1,1,1)), dtype=float)
    try:    ### This is a very dangerous conditional statement to put
        binmeth_fiter = np.array(binmeth_whole)[:,[2,4,6]]  ## filter for columns from thw binmeth_whole
        binmeth_fiter[binmeth_fiter == -1] = np.nan  #### Removing the reads which do not have some information
        binmeth_fiter = binmeth_fiter[~np.isnan(binmeth_fiter).any(axis=1)]
        kmeans = KMeans(n_clusters=8, init=init_cls,n_init = 1).fit(binmeth_fiter)
        type_cls = np.unique(kmeans.labels_, return_counts=True)
        return(countTypeFreqs(type_cls))
    except:  ### Just skipped the value error, where all zero rows and columns in binmeth_fiter
        type_counts = [0,0,0,0,0,0,0,0]
        type_freqs = [0,0,0,0,0,0,0,0]
        return(type_counts,type_freqs)

#### ==========================================
##      MHL calculation from here
##  adapted from Guo et al, 2017
#### ==========================================


def getCs_bins_alignment(bins_alignment, strand):
    # re_string = ['C','T','.']
    re_string = getStrandContext(strand)
    num_seq = len(bins_alignment)
    refseq = bins_alignment[0]
    cs_hap_bins = []
    for seqind in range(1, num_seq):
        eseq = str(bins_alignment[seqind].seq)
        eseq = eseq.replace(re_string[2], "")
        eseq = eseq.replace(re_string[0], "1")
        eseq = eseq.replace(re_string[1], "0")
        if len(eseq) > 0:
            cs_hap_bins.append(eseq)
    return(np.array(cs_hap_bins))

def getCs_seq_record(seqrecord, strand):
    re_string = getStrandContext(strand)
    eseq = str(seqrecord.seq)
    eseq = eseq.replace(re_string[2], "")
    eseq = eseq.replace(re_string[0], "1")
    eseq = eseq.replace(re_string[1], "0")
    ## Sometimes this is happening,  '10T00000' when there is a SNP?
    return(eseq)

def get_hap_comb_single(ecs_hap_bins):
    ncols = len(ecs_hap_bins)
    hap_comb_single = np.zeros(0, dtype=ecs_hap_bins.dtype)
    for i in range(ncols):
        for j in range(i, ncols):
            hap_comb_single = np.append(hap_comb_single, ecs_hap_bins[i:j+1])
    return(hap_comb_single)

def get_hap_comb(cs_hap_bins):
    hap_comb = np.zeros(0, dtype=cs_hap_bins.dtype)
    for i in range(len(cs_hap_bins)):
        hap_comb = np.append(hap_comb, get_hap_comb_single(cs_hap_bins[i]))
    return(hap_comb)

def calculate_mhl(cs_hap_bins):
    ## this is adapted from Guo et al. 2017 Figure 2.
    if len(cs_hap_bins) <= 1:
        return((np.nan, 0, 0))
    nplen = np.vectorize(len)
    hap_comb = get_hap_comb(cs_hap_bins)
    mhl_value = np.zeros(0, dtype=float)
    ncols = np.max(nplen(hap_comb))
    for i in range(ncols):
        input_tmp = np.where(nplen(hap_comb) == i + 1)[0]
        nmeth = len( np.where( hap_comb[input_tmp] == '1' * (i + 1) )[0]  )
        if len(input_tmp) >= len(cs_hap_bins)/1.3:
            mhl_value = np.append(mhl_value, (i + 1) * nmeth / float(len(input_tmp)))
    mhl_value_final = np.sum(mhl_value) / sum(range(ncols + 1))
    return((mhl_value_final, len(cs_hap_bins), ncols))

def get_mhl_entire_bed(inBam, meths, tair10, required_bed, outstat = ''):
    window_size = required_bed[3]
    read_length_thres = window_size/2
    bin_start = required_bed[1] - (required_bed[1] % window_size)
    estimated_bins = list(range(bin_start, required_bed[2], window_size))
    progress_bins = 0
    #### iter meths in the required bed.
    meths_bins = meths.iter_bed_windows([required_bed[0], bin_start, required_bed[2]], window_size)
    for bins, temp_meth_bins in zip(estimated_bins, meths_bins) :        ## sliding windows with window_size
        progress_bins += 1
        cs_hap_bins = np.zeros(0, dtype="string")
        bin_bed = [required_bed[0], bins, bins + window_size - 1]
        for binread in inBam.fetch(str(bin_bed[0]), bin_bed[1], bin_bed[2]):
            r_strand = getStrandPE(binread)
            if filterRead(binread): ## Filtering the dupplicated reads.
                continue
            intersect_bed = [bin_bed[0], binread.reference_start, binread.reference_end]
            rseq_record = getSeqRecord(binread, tair10, intersect_bed, [intersect_bed[1], intersect_bed[2]])
            cs_bin = getCs_seq_record(rseq_record, r_strand)
            if len(cs_bin) > 0:
                cs_hap_bins = np.append(cs_hap_bins, cs_bin)
        mhl_value, no_cs_hap_bins, ncols = calculate_mhl(cs_hap_bins)
        wma_win = meth5py.MethylationSummaryStats(meths, temp_meth_bins[1], 1)
        if wma_win == 0 or np.isnan(wma_win):
            frac_mhl = mhl_value
        else:
            frac_mhl = mhl_value
        if outstat == '':
            print(("%s,%s,%s,%s,%s,%s" % (bin_bed[0], str(bin_bed[1]), str(bin_bed[2]), frac_mhl, wma_win, no_cs_hap_bins)))
        else:
            outstat.write("%s,%s,%s,%s,%s,%s\n" % (bin_bed[0], str(bin_bed[1]), str(bin_bed[2]), frac_mhl, wma_win, no_cs_hap_bins))
        if progress_bins % 1000 == 0:
            log.info("ProgressMeter - %s windows in analysed, %s total" % (progress_bins, len(estimated_bins)))
    return(0)

def potato_mhl_calc(args):
    log.info("loading the input files!")
    inBam = pysam.AlignmentFile(args['inFile'], "rb")
    if args['inhdf5'] != '':
        meths = meth5py.load_hdf5_methylation_file(args['inhdf5'])
    else:
        meths = ''
    (chrs, chrslen, binLen) = getChrs(inBam)
    tair10 = Fasta(args['fastaFile'])
    if args['window_size'] is None:
        window_size = binLen
    else:
        window_size = args['window_size']
    log.info("done!")
    if args['outFile'] != "STDOUT":
        outstat = open(args['outFile'], 'w')
        outstat.write("chr,start,end,mhl_stat,wma_win,counts\n")
    else:
        outstat = ''
        print("chr,start,end,mhl_stat,wma_win,counts")
    if args['reqRegion'] != '0,0,0':
        required_region = args['reqRegion'].split(',')
        required_bed = [required_region[0], int(required_region[1]), int(required_region[2]), window_size]
        log.info("analysing region %s:%s-%s !" % (required_bed[0], required_bed[1], required_bed[2]))
        get_mhl_entire_bed(inBam, meths, tair10, required_bed, outstat)
        log.info("finished!")
        return(0)
    for cid, clen in zip(chrs, chrslen):     ## chromosome wise
        log.info("analysing chromosome: %s" % cid)
        required_bed = [cid, 1, clen, window_size]   ### 0 is to not print reads
        get_mhl_entire_bed(inBam, meths, tair10, required_bed, outstat)
    log.info("finished!")
    return(0)

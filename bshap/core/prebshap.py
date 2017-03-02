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
    return (reqchrs, reqchrslen, int(np.nanmean(binLen)))

def get_reverse_complement(seq):
    old_chars = "ACGT"
    replace_chars = "TGCA"
    tab = string.maketrans(old_chars,replace_chars)
    return seq.translate(tab)[::-1]

def findIntersectbed(bin_bed, map_bed):
    # bin_bed = [290605, 290688]
    # map_bed = [290626,290711]
    bin_bed_arr = np.array(range(bin_bed[0], bin_bed[1]))
    map_bed_arr = np.array(range(map_bed[0], map_bed[1]))
    intersect_bed_arr = np.sort(np.intersect1d(bin_bed_arr, map_bed_arr))
    intersect_bed = [intersect_bed_arr[0], intersect_bed_arr[-1]+1]
    return intersect_bed

def parseContext(tair10, cid, pos, strand):
    try:
        if strand == '0':
            dnastring = tair10[cid][pos:pos+3].seq.encode('ascii').upper()
        elif strand == '1': ## make sure you can have to identify strand here
            dnastring = tair10[cid][pos-2:pos+1].seq.encode('ascii').upper()  ## Changed the position, different from forward
            dnastring = get_reverse_complement(dnastring)
    except:
        return ("CNN", "CN")
    if dnastring[1].upper() == 'G':
        dna_context = ["CG",0]
    elif dnastring[2].upper() == 'G':
        dna_context = ["CHG",1]
    elif dnastring:
        dna_context = ["CHH",2]
    return (dnastring, dna_context)

def decodeFlag(flag):
    bintag = str(int(bin(flag)[2:]))[::-1]
    notag = (12 - len(bintag)) * '0'
    return bintag + notag

def getHighlightedSeqs(refseq, rseq, strand):
    if strand == '0': ## forward, let C's stay
        re_string = ['C','T','.']
    elif strand == '1':  ## G's should stay
        re_string = ['G', 'A',',']
    dot_rseq = rseq
    for i,c in enumerate(refseq):
        if c.upper() != re_string[0]:
            dot_rseq = dot_rseq[:i] + re_string[2] + dot_rseq[i+1:]
    return dot_rseq

def getMatchedSeq(binread, tair10, bin_bed):
    intersect_bed = findIntersectbed(bin_bed[1:3], [binread.reference_start, binread.reference_end])
    bin_positions = np.where(np.in1d(binread.get_reference_positions(), range(intersect_bed[0], intersect_bed[1])))[0]
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
    return refseq, rseq

def getSeqRecord(binread, tair10, bin_bed, intersect_bed):
    ## bin_bed = ['Chr1', start, end, binLen, pointPos]
    ## intersect_bed = [s, e]
    char_add = 'N'
    strand = decodeFlag(binread.flag)[4]
    refseq, rseq = getMatchedSeq(binread, tair10, bin_bed)
    dot_rseq = getHighlightedSeqs(refseq, rseq, strand)
    fin_rseq = char_add * (intersect_bed[0] - bin_bed[1]) + dot_rseq + char_add * (bin_bed[2] - intersect_bed[1])
    return SeqRecord(Seq(fin_rseq, generic_dna), id = binread.query_name.split(' ')[0])

def getMethRead(tair10, binread):
    strand = decodeFlag(binread.flag)[4]
    error_rate = 0.01
    if strand == '0': ## Forward read, count the number of C's
        context = "C"
        #reqstr = "T"   # if this is present the base is non methylated
    elif strand == '1': ## reverse
        context = "G"
        #reqstr = "A"   # if this is present the base is non methylated
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
    nc_thres = 0
    rmeth = [float(mCs[i])/tCs[i] if tCs[i] > nc_thres else -1 for i in range(3)]
    tmeth = float(sum(mCs))/sum(tCs) if sum(tCs) > nc_thres else -1
    permeths = [tmeth, rmeth[0], rmeth[1], rmeth[2]]
    return permeths
    # We need to differentiate between the reads which are small and has no methylation

def getMethWind(inBam, tair10, required_bed, meths):
    # required_bed = ['Chr1', start, end, binLen]
    binLen = required_bed[3]
    read_length_thres = binLen/2
    bin_start = required_bed[1] - (required_bed[1] % binLen)
    estimated_bins = range(bin_start, required_bed[2], binLen)
    dt = h5py.special_dtype(vlen=np.dtype('float16'))
    reqmeths = meths.create_dataset(required_bed[0],compression="gzip", shape = (len(estimated_bins),4,), dtype=dt)
    meths[required_bed[0]].attrs['positions'] = estimated_bins
    progress_bins = 0
    window_alignment = []  ## multiple segment alignment
    for bins in estimated_bins:        ## sliding windows with binLen
        binmeth = []
        progress_bins += 1
        dot_refseq = tair10[required_bed[0]][bins:bins + binLen].seq.encode('ascii')
        dot_refseq = re.sub('A|T', '.', dot_refseq.upper())
        bins_alignment = [SeqRecord(Seq(dot_refseq, generic_dna), id = required_bed[0] + ':' + str(bins) + '-' + str(bins + binLen))]
        for binread in inBam.fetch(required_bed[0], bins, bins + binLen):
            rflag = decodeFlag(binread.flag)
            bin_bed = [required_bed[0], bins, bins + binLen]
            if len(binread.get_aligned_pairs()) == 0 and rflag[10] != '0' and rflag[8] != '0': ## removing the duplicates
                continue        ## Skip the read if the alignment is not good ;)
            intersect_bed = findIntersectbed(bin_bed[1:3], [binread.reference_start, binread.reference_end])
            intersect_len = len(range(intersect_bed[0], intersect_bed[1]))
            rseq_record = getSeqRecord(binread, tair10, bin_bed, intersect_bed) ## No need to check for the overlap to print the alignment
            if intersect_len > read_length_thres:
                permeths = getMethRead(tair10, binread) ## taking the first and last
                binmeth.append(permeths)
            bins_alignment.append(rseq_record)
        if progress_bins % 1000 == 0:
            log.info("ProgressMeter - %s windows in analysed, %s total" % (progress_bins, estimated_bins))
        reqmeths[progress_bins-1] = np.array(binmeth).T
        window_alignment.append(MultipleSeqAlignment(bins_alignment))
    return window_alignment

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
            required_bed = [cid, 0, clen, binLen]   ### 0 is to not print reads
            window_alignment = getMethWind(inBam, tair10, required_bed, meths)
    else:
        required_bed = interesting_region.split(',')
        if len(required_bed) == 3:
            required_bed = [required_bed[0], int(required_bed[1]), int(required_bed[2]), binLen, 0]
        else:
            die("error in the input bed position: chromosome,start,end")
        log.info("analysing region %s:%s-%s !" % (required_bed[0], required_bed[1], required_bed[2]))
        window_alignment = getMethWind(inBam, tair10, required_bed, meths)
        AlignIO.write(window_alignment, 'meths.' + outFile + '.aln', "clustal")
    meths.close()
    log.info("finished!")

def getMethReadsgivenPosition(bamFile):
    return 0

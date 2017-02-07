## Main module for pyBsHap

import pysam
from pyfaidx import Fasta
import logging
import numpy as np
import re
import json
import string

log = logging.getLogger(__name__)

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
        dna_context = "CG"
    elif dnastring[2].upper() == 'G':
        dna_context = "CHG"
    elif dnastring:
        dna_context = "CHH"
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
    fin_refseq = re.sub('A|T', re_string[2], refseq.upper()).lower()
    fin_refseq = re.sub(re_string[0].lower(), re_string[0], fin_refseq)
    fin_rseq = rseq
    for i,c in enumerate(refseq):
        if c.upper() != re_string[0]:
            fin_rseq = fin_rseq[:i] + re_string[2] + fin_rseq[i+1:]
    return fin_refseq, fin_rseq

def printInterestingRegion(tair10, binread, rmeth, pointPos=0):
    ## pointPos = 290598
    ## map_pos = [cid, ref_pos_start, ref_pos_start_end]
    map_pos = [binread.reference_name, binread.reference_start, binread.reference_end]
    strand = decodeFlag(binread.flag)[4]
    refseq = tair10[map_pos[0]][map_pos[1]:map_pos[2]].seq.encode('ascii')
    rseq = binread.seq
    dot_refseq, dot_rseq = getHighlightedSeqs(refseq, rseq, strand)
    print("start:end:meth -- %s: %s - %s : %s" % (map_pos[0], map_pos[1], map_pos[2], rmeth))
    if pointPos < map_pos[2] and pointPos > map_pos[1]:
        point_len = pointPos - map_pos[1]
        print("%s" % '.' * point_len + `pointPos`)
    print("%s" % dot_refseq)
    print("%s" % dot_rseq)

def getMethRead(reqcontext, tair10, binread):
    map_pos = [binread.reference_name, binread.reference_start, binread.reference_end]
    refseq = tair10[map_pos[0]][map_pos[1]:map_pos[2]].seq.encode('ascii')
    rseq = binread.seq
    strand = decodeFlag(binread.flag)[4]
    error_rate = 0.01
    if strand == '0': ## Forward read, count the number of C's
        context = "C"
        #reqstr = "T"   # if this is present the base is non methylated
    elif strand == '1': ## reverse
        context = "G"
        #reqstr = "A"   # if this is present the base is non methylated
    if len(refseq) != len(rseq):
        return "check the strings!"
    mc, mt, read_length = 0, 0, 0
    for i,c in enumerate(refseq):
        read_length = read_length + 1
        if c.upper() == context:
            (dnastring, methcontext) = parseContext(tair10, map_pos[0], map_pos[1] + i, strand)
            checkcontext = reqcontext == methcontext
            if reqcontext == 'CN':
                checkcontext = True
            if checkcontext:
                mt += 1
                if rseq[i].upper() == context:
                    mc += 1
    #log.debug("%s:%s:%s:%s" % (ref_pos,mc, mt, read_length))
    if mt >= 1:
        return mc, mt, float(mc)/mt
    else:
        return mc, mt, -1   ### We need to differentiate between the reads which are small and has no methylation

def getMethWind(inBam, tair10, bin_bed, reqcontext = 'CN'):
    # bin_bed = ['Chr1', start, end, binLen, pointPos]
    binLen = bin_bed[3]
    read_length_thres = binLen/2
    meths = {}
    for bins in range(bin_bed[1], bin_bed[2], binLen):        ## sliding windows with binLen
        binmeth = np.zeros(0)
        binmc = np.zeros(0)
        binmt = np.zeros(0)
        for binread in inBam.fetch(bin_bed[0], bins, bins + binLen):
            rflag = decodeFlag(binread.flag)
            intersect_bed = findIntersectbed([bins, bins + binLen], [binread.reference_start, binread.reference_end])
            intersect_len = len(range(intersect_bed[0], intersect_bed[1]))
            if intersect_len > read_length_thres and rflag[10] == '0' and rflag[8] == '0': ## removing the duplicates
                ## The start of the read is binread.reference_start
                mc, mt, rmeth = getMethRead(reqcontext, tair10, binread) ## taking the first and last
                if not np.isnan(rmeth):
                    binmeth = np.append(binmeth, rmeth)
                    binmc = np.append(binmc, mc)
                    binmt = np.append(binmt, mt)
                    if bin_bed[4] != 0:
                        printInterestingRegion(tair10, binread, rmeth, bin_bed[4])
        meths[bins+1] = binmeth.tolist()
    return meths

def getMethGenome(bamFile, fastaFile, outFile, reqcontext = "CN", interesting_region='0,0,0'):
    ## input a bam file, bin length
    ## make sure the binLen is less than the read length
    inBam = pysam.AlignmentFile(bamFile, "rb")
    (chrs, chrslen, binLen) = getChrs(inBam)
    tair10 = Fasta(fastaFile)
    if interesting_region == '0,0,0,0':
        meths = {}
        meths["chrs"] = chrs
        meths["binlen"] = binLen
        meths["chrslen"] = chrslen
        for cid, clen in zip(chrs, chrslen):     ## chromosome wise
            log.info("analysing chromosome: %s" % cid)
            bin_bed = [cid, 0, clen, binLen, 0]   ### 0 is to not print reads
            meths[cid] = getMethWind(inBam, tair10, bin_bed, reqcontext=reqcontext)
        with open(outFile, 'wb') as fp:
            fp.write(json.dumps(meths))
    else:
        meths = {}
        bin_bed = interesting_region.split(',')
        if len(bin_bed) == 3:
            bin_bed = [bin_bed[0], int(bin_bed[1]), int(bin_bed[2]), binLen, 0]
        else:
            bin_bed = [bin_bed[0], int(bin_bed[1]), int(bin_bed[2]), binLen, int(bin_bed[3])]
        meths[bin_bed[0]] = getMethWind(inBam, tair10, bin_bed, reqcontext)
        with open(outFile, 'wb') as fp:
            fp.write(json.dumps(meths))

import pysam
from pyfaidx import Fasta
import logging
import numpy as np
import re
import json

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

def parseContext(tair10, cid, pos):
    try:
        dnastring = tair10[cid][pos:pos+3].seq.encode('ascii')
    except:
        return None
    if dnastring[1].upper() == 'G':
        return "CG"
    elif dnastring[2].upper() == 'G':
        return "CHG"
    else:
        return "CHH"

def decodeFlag(flag):
    bintag = str(int(bin(flag)[2:]))[::-1]
    notag = (12 - len(bintag)) * '0'
    return bintag + notag

def getMethRead(reqcontext, tair10, cid, bins, refseq, rseq, strand):
    error_rate = 0.01
    if strand == '0': ## Forward read, count the number of C's
        context = "C"
        reqstr = "T"
    elif strand == '1': ## reverse
        context = "G"
        reqstr = "A"
    if len(refseq) != len(rseq):
        return "check the strings!"
    mc, mt = 0, 0
    for i,c in enumerate(refseq):
        methcontext = parseContext(tair10, cid, bins + i)
        checkcontext = reqcontext == methcontext
        if reqcontext == 'CN':
            checkcontext = True
        if c.upper() == context and checkcontext:
            if rseq[i].upper() == context:
                mc += 1
                mt += 1
            elif rseq[i].upper() == reqstr:
                mt += 1
    if mt > 5:
        return mc, mt, float(mc)/mt
    else:
        return mc, mt, np.nan

def getMethWind(bamFile, fastaFile, reqcontext, outFile):
    ## input a bam file, bin length
    ## make sure the binLen is less than the read length
    seqthres = 20 ## minimum overlap within the window to get the required
    inBam = pysam.AlignmentFile(bamFile, "rb")
    (chrs, chrslen, binLen) = getChrs(inBam)
    tair10 = Fasta(fastaFile)
    meths = {}
    meths["chrs"] = chrs
    meths["binlen"] = binLen
    meths["chrslen"] = chrslen
    for cid, clen in zip(chrs, chrslen):     ## chromosome wise
        log.info("analysing chromosome: %s" % cid)
        meths[cid] = {}
        for bins in range(0, clen, binLen):        ## sliding windows with binLen
            binmeth = np.zeros(0)
            binmc = np.zeros(0)
            binmt = np.zeros(0)
            for binread in inBam.fetch(cid, bins, bins + binLen):
                rflag = decodeFlag(binread.flag)
                rind = np.where(np.in1d(np.array(binread.get_reference_positions()), np.array(range(bins,bins + binLen))))[0]
                refseq = tair10[cid][binread.reference_start:binread.reference_end].seq.encode('ascii')
                if len(rind) > seqthres and rflag[10] == '0' and rflag[8] == '0': ## removing the duplicates
                    mc, mt, rmeth = getMethRead(reqcontext, tair10, cid, bins, refseq[rind[0]:rind[-1]], binread.seq[rind[0]:rind[-1]], rflag[4]) ## taking the first and last
                    if not np.isnan(rmeth):
                        binmeth = np.append(binmeth, rmeth)
                        binmc = np.append(binmc, mc)
                        binmt = np.append(binmt, mt)
            meths[cid][bins+1] = binmeth.tolist()
    with open(outFile, 'wb') as fp:
        fp.write(json.dumps(meths))

#def getMethReads(bamFile, fastaFile, outFile):

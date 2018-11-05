import pysam
from pyfaidx import Fasta
import logging
import numpy as np
import os.path, sys
from subprocess import Popen, PIPE
import shlex

from . import prebshap

log = logging.getLogger(__name__)

def die(msg):
  sys.stderr.write('Error: ' + msg + '\n')
  sys.exit(1)

def md_counter(fin_rseq, re_string):
    md_tag = ''
    check_str = fin_rseq
    while len(check_str) > 0:
        strip_num = len(check_str) - len(check_str.lstrip(re_string[2]))
        if strip_num > 0:
            md_tag = md_tag + str(strip_num)
            check_str = check_str.lstrip(re_string[2])
        else:
            md_tag = md_tag + re_string[0]
            check_str = check_str[1:]
    return md_tag

def modifyMDtag(inBam, tair10, binread, outBam):
    ## Modifying MD tag in reads for jbrowse viewing
    ##
    ## Need to define an outBam before calling this command
    ##with pysam.AlignmentFile(outBam_Name, "wb", header = inBam.header) as outBam:
    ##
    strand = prebshap.decodeFlag(binread.flag)[4]
    oread = pysam.AlignedSegment()
    oread = binread
    tags = [i for i in binread.get_tags() if i[0] != 'MD']
    req_bin = [binread.reference_start, binread.reference_end]
    fin_rseq = prebshap.getSeqRecord(binread, tair10, [binread.reference_name, req_bin[0], req_bin[1]], req_bin, defreturn=1)
    re_string = prebshap.getStrandContext(strand)
    md_tag = md_counter(fin_rseq, re_string)
    tags.append(('MD', md_tag))
    oread.tags = tags
    outBam.write(oread)

def potatoskin_modify_mdtag_bam(bamFile, fastaFile, outFile):
    # Out Bam file
    log.info("loading the input files!")
    inBam = pysam.AlignmentFile(bamFile, "rb")
    (chrs, chrslen, binLen) = prebshap.getChrs(inBam)
    tair10 = Fasta(fastaFile)
    log.info("finished!")
    log.info("writing file into AlignmentFile, %s!" % outFile)
    outBam = pysam.AlignmentFile(outFile, "wb", header = inBam.header)
    for cid, clen in zip(chrs, chrslen):     ## chromosome wise
        log.info("analysing chromosome: %s" % cid)
        for binread in inBam.fetch(str(cid), 0, clen):
            if prebshap.filterRead(binread):
                continue
            modifyMDtag(inBam, tair10, binread, outBam)
        log.info("finished!")
    log.info("indexing output bam file!")
    Popen(shlex.split("samtools index " + outFile), stdout=PIPE)
    log.info("finished!")

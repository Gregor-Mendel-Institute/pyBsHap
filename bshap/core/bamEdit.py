import pysam
from pyfaidx import Fasta
import logging
import numpy as np
import re
import string
import os.path, sys
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from subprocess import Popen, PIPE
import shlex

import prebshap

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

def writeBam(bamFile, fastaFile, outFile, interesting_region='0,0,0'):
    # Out Bam file
    log.info("loading the input files!")
    inBam = pysam.AlignmentFile(bamFile, "rb")
    (chrs, chrslen, binLen) = prebshap.getChrs(inBam)
    tair10 = Fasta(fastaFile)
    log.info("finished!")
    outBam_file = outFile + '_processed_reads_no_clonal_modifiedMD_filtered.bam'
    log.info("writing file into AlignmentFile, %s!" % outBam_file)
    outBam = pysam.AlignmentFile(outBam_file, "wb", header = inBam.header)
    if interesting_region == '0,0,0':
        for cid, clen in zip(chrs, chrslen):     ## chromosome wise
            log.info("analysing chromosome: %s" % cid)
            for binread in inBam.fetch(str(cid), 0, clen):
                if prebshap.filterRead(binread):
                    continue
                modifyMDtag(inBam, tair10, binread, outBam)
            log.info("finished!")
    else:
        required_region = interesting_region.split(',')
        required_bed = [required_region[0], int(required_region[1]), int(required_region[2]), binLen, 1]
        log.info("analysing region %s:%s-%s !" % (required_bed[0], required_bed[1], required_bed[2]))
        for binread in inBam.fetch(required_bed[0], required_bed[1], required_bed[2]):
            if prebshap.filterRead(binread):
                continue
            modifyMDtag(inBam, tair10, binread, outBam)
        log.info("finished!")
    log.info("indexing output bam file!")
    Popen(shlex.split("samtools index " + outBam_file), stdout=PIPE)
    log.info("finished!")

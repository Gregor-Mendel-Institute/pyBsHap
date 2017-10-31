"""
Pakaged functions to get the methyaltion levels using bed tools
"""
import os.path
import numpy as np
from subprocess import Popen, PIPE
import logging
import itertools

log = logging.getLogger(__name__)
chrs = ['Chr1','Chr2','Chr3','Chr4','Chr5']
golden_chrlen = [30427671, 19698289, 23459830, 18585056, 26975502]
entire_chrslen = [34964571, 22037565, 25499034, 20862711, 31270811]

def windows(seqlength, window_size, overlap):
    if overlap >= window_size:
        raise NotImplementedError
    if overlap > 0:
        for x in range(1, seqlength, overlap):
            yield([x, x + window_size - 1])
    else:
        for x in range(1, seqlength, window_size):
            yield([x, x + window_size - 1])

def generate_window_file(window_size, out_windows, overlap):
    if os.path.isfile(out_windows):
        log.info("utilizing window file: %s" % out_windows)
        return 0
    log.info("generating a window file: %s" % out_windows)
    outWindow = open(out_windows, 'w')
    for echr, echrlen in zip(chrs, golden_chrlen):
        echr_windows = windows(echrlen, window_size, overlap)
        for ewind in echr_windows:
            outWindow.write('%s\t%s\t%s\n' %(echr, ewind[0], ewind[1]))
    outWindow.close()
    log.info("done!")
    return 0

def get_genomewide_methylation_WeightedMean(bedtoolPath, bedFile, outFile, window_size, overlap):
    outBedGraph = open(outFile, "w")
    window_file = "tair10." + str(window_size) + "bp_windowsize." + str(overlap) + "bp_overlap.windows.txt"
    generate_window_file(window_size, window_file, overlap)
    bedtools_command = 'bedtools map -a ' + window_file + ' -b ' + bedFile + ' -o sum,sum -c 7,8'
    if bedtoolPath is not None:
        bedtools_command = bedtoolPath + '/' + bedtools_command
    awk_command = '| awk \'$5 > 0 {print $1 "\t" $2 "\t" $3 "\t" $4/$5}\''
    full_command = bedtools_command + awk_command
    log.info('running bedtools!')
    convertcsv = Popen(full_command, shell=True, stdout = outBedGraph)
    convertcsv.wait()
    outBedGraph.close()
    log.info('done!')

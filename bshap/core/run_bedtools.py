"""
Pakaged functions to get the methyaltion levels using bed tools
"""
import os.path
import numpy as np
from subprocess import Popen, PIPE
import logging

log = logging.getLogger(__name__)
chrs = ['Chr1','Chr2','Chr3','Chr4','Chr5']
golden_chrlen = [30427671, 19698289, 23459830, 18585056, 26975502]
entire_chrslen = [34964571, 22037565, 25499034, 20862711, 31270811]

def generate_window_file(window_size, out_windows):
    if os.path.isfile(out_windows):
        return 0
    outWindow = open(out_windows, 'w')
    for echr, echrlen in zip(chrs, golden_chrlen):
        echr_windows = range(1, echrlen, window_size)
        for ewind in echr_windows:
            outWindow.write('%s\t%s\t%s\n' %(echr, ewind, ewind + window_size - 1))
    outWindow.close()
    return 0

def get_genomewide_methylation_WeightedMean(bedtoolPath, bedFile, outFile, window_size=100):
    outBedGraph = open(outFile, "w")
    window_file = "temp." + str(window_size) + "bp.windows.txt"
    generate_window_file(window_size, window_file)
    bedtools_command = 'bedtools map -a ' + window_file + ' -b ' + bedFile + ' -o sum,sum -c 7,8'
    if bedtoolPath != '':
        bedtools_command = bedtoolPath + '/' + bedtools_command
    awk_command = '| awk \'$5 > 0 {print $1 "\t" $2 "\t" $3 "\t" $4/$5}\''
    full_command = bedtools_command + awk_command
    log.info('running bedtools!')
    convertcsv = Popen(full_command, shell=True, stdout = outBedGraph)
    convertcsv.wait()
    outBedGraph.close()
    log.info('done!')

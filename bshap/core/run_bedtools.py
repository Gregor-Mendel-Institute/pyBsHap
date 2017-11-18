"""
Pakaged functions to get the methyaltion levels using bed tools
"""
import os.path
import numpy as np
from subprocess import Popen, PIPE
import logging
import itertools
import pybedtools

log = logging.getLogger(__name__)
chrs = ['Chr1','Chr2','Chr3','Chr4','Chr5']
golden_chrlen = [30427671, 19698289, 23459830, 18585056, 26975502]
entire_chrslen = [34964571, 22037565, 25499034, 20862711, 31270811]

def windows(seqlength, window_size, overlap):
    if overlap >= window_size:
        raise(NotImplementedError)
    if overlap > 0:
        for x in range(1, seqlength, overlap):
            yield([x, x + window_size - 1])
    else:
        for x in range(1, seqlength, window_size):
            yield([x, x + window_size - 1])

def generate_window_file(window_size, overlap, out_windows):
    if os.path.isfile(out_windows):
        log.info("utilizing window file: %s" % out_windows)
        return(0)
    log.info("generating a window file: %s" % out_windows)
    outWindow = open(out_windows, 'w')
    for echr, echrlen in zip(chrs, golden_chrlen):
        echr_windows = windows(echrlen, window_size, overlap)
        for ewind in echr_windows:
            outWindow.write('%s\t%s\t%s\n' %(echr, ewind[0], ewind[1]))
    outWindow.close()
    log.info("done!")
    return(0)

def get_weighted_average(feature):
    """
    Parsing the output of bedtools to get weighted average
    """
    if feature[3] != '.':
        if int(feature[3]) > 0:
            num_c = np.sum(np.array(feature[4].split(","), dtype=int) * np.array(feature[5].split(","), dtype=float))
            return(num_c / int(feature[3]))
    return(None)

def MethylationSummaryStats(window_file, bedFile, bedtoolsPath, category):
    ## In a bedfile you have these columns
    # Chr1	22	23	CCC	0	+	1	1
    # Chr1	23	24	CCT	0	+	1	1
    # Chr1	24	25	CTA	0	+	1	1
    # Chr1	29	30	CCT	0	+	0	1
    # Chr1	30	31	CTC	0	+	1	1
    # Chr1	32	33	CTG	0	+	1	1
    # Chr1	34	35	CAG	1	-	9	18
    bedtools_command = 'bedtools map -a ' + window_file + ' -b ' + bedFile
    if bedtoolsPath is not None:
        bedtools_command = bedtoolsPath + '/' + bedtools_command
    skip_na_lines = ' | awk \'$4 != "." {print $0}\''
    if category == 1:   # weighted mean
        bedtools_command = bedtools_command + ' -o sum,collapse,collapse -c 8,5,7'
        return(bedtools_command)
    elif category == 2:  # fraction of methylated positions
        bedtools_command = bedtools_command + ' -o mean -c 5'
        return(bedtools_command + skip_na_lines)
    elif category == 3:  # absolute means
        awk_command = 'awk \'$8 > 0 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $7/$8 "\t" $6}\' ' + bedFile
        bedtools_command =  'bedtools map -a ' + window_file + ' -b stdin -o mean -c 5'
        if bedtoolsPath is not None:
            bedtools_command = bedtoolsPath + '/' + bedtools_command
        return(awk_command + ' | ' + bedtools_command + skip_na_lines)
    else:
        raise(NotImplementedError)

def get_genomewide_methylation_WeightedMean(args):
    category = args['category']
    outBedGraph = open(args['outFile'], "w")
    window_file = "tair10." + str(args['window_size']) + "bp_windowsize." + str(args['overlap']) + "bp_overlap.windows.txt"
    generate_window_file(args['window_size'], args['overlap'], window_file)
    full_command = MethylationSummaryStats(window_file, args['inFile'], args['bedtoolsPath'], args['category'])
    log.info("make sure the bedtools version is > v2.26.0")
    log.info('running bedtools!')
    if category != 1:
        convertcsv = Popen(full_command, shell=True, stdout = outBedGraph)
        convertcsv.wait()
        outBedGraph.close()
        log.info("done!")
    elif category == 1:
        convertcsv = Popen(full_command, shell=True, stdout = PIPE)
        for each_line in iter(convertcsv.stdout.readline, ''):
            each_bed = each_line.rstrip().split("\t")
            score = get_weighted_average(each_bed)
            if score is not None:
                outBedGraph.write("%s\t%s\t%s\t%s\n" % (each_bed[0], each_bed[1], each_bed[2], score))
        outBedGraph.close()
        log.info('done!')

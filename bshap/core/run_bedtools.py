"""
Pakaged functions to get the methyaltion levels using bed tools
"""
import os.path
import numpy as np
from subprocess import Popen, PIPE
import logging
import pandas as pd

log = logging.getLogger(__name__)
chrs = ['Chr1','Chr2','Chr3','Chr4','Chr5']
golden_chrlen = [30427671, 19698289, 23459830, 18585056, 26975502]
entire_chrslen = [34964571, 22037565, 25499034, 20862711, 31270811]

def identify_positions_given_names(in_file, araport11_file):
    if araport11_file is None:
        raise NameError("please give araport file to get positions")
    araport11 = pd.read_csv(araport11_file, header = None, sep = "\t")
    bed_names = np.array(pd.read_csv(in_file, header = None)[0])
    req_bed_df = araport11.loc[araport11[3].isin(bed_names),]
    return(req_bed_df)

def get_filter_bed_ix(bed_file, input_bed, just_names=True, araport11_file=None):
    ## Here input_bed is a pandas dataframe
    import pybedtools as pybed
    if not isinstance(input_bed, pd.DataFrame):
        raise(NotImplementedError)
    ## In case of just_names the the bed_file contains only names
    if not os.path.isfile(bed_file):
        raise NameError("file is not present, check file name")
    if just_names:
        req_bed_df = identify_positions_given_names(bed_file, araport11_file)
        req_bed = pybed.BedTool.from_dataframe(req_bed_df.loc[:,[0,1,2]])
    else:
        req_bed = pybed.BedTool(bed_file)
    inBed = pybed.BedTool.from_dataframe(input_bed.loc[:,[0,1,2]])
    ## Just taking first three columns for bedtools
    unionBed = inBed.intersect(req_bed, wa=True).to_dataframe() ## wa is to return the entire bed.
    total_cols = np.array(input_bed.loc[:,0] + "," + input_bed.loc[:,1].map(str) + "," +  input_bed.loc[:,2].map(str), dtype="str")
    unionBed_cols = np.array(unionBed["chrom"] + "," + unionBed["start"].map(str) + "," +  unionBed["end"].map(str), dtype="str")
    return(np.where(np.in1d(total_cols, unionBed_cols))[0])

def get_filter_pos_echr(bed_file, chrid, common_positions, just_names = True, araport11_file=None):
    filter_inds = []
    if just_names:
        req_name_pos = identify_positions_given_names(bed_file, araport11_file)
        req_name_pos = req_name_pos.loc[req_name_pos[0] == chrid]
    else:
        req_name_pos = pd.read_table(bed_file, header = None)
        req_name_pos = req_name_pos.loc[req_name_pos[0] == chrid]
    for e_ind, e_pos in req_name_pos.iterrows():
        pos_inds = np.searchsorted(common_positions,[e_pos[1],e_pos[2]], side='right')
        filter_inds.extend(range(pos_inds[0],pos_inds[1]))
        #color_scatter[pos_inds[0]:pos_inds[1]] = replace_to
    return(np.array(filter_inds))

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
    if len(feature) < 6:
        raise(NotImplementedError, "the bed file provided should have atleast 3 columns")
    if feature[-3] != '.':
        if int(feature[-3]) > 0:
            num_c = np.sum(np.array(feature[-2].split(","), dtype=int) * np.array(feature[-2].split(","), dtype=float))
            return(num_c / int(feature[-3]))
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
    genome_out = open("tair10.genome.txt", 'w')
    for echr, echrlen in zip(chrs, golden_chrlen):
        genome_out.write("%s\t%s\n" % (echr, echrlen))
    genome_out.close()
    bedtools_command = 'bedtools map -g tair10.genome.txt -a ' + window_file + ' -b ' + bedFile
    if bedtoolsPath is not None:
        bedtools_command = bedtoolsPath + '/' + bedtools_command
    skip_na_lines = ' | awk \'$4 != "bedtools_command." {print $0}\''
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

def get_genomewide_methylation_average(args):
    category = args['category']
    outBedGraph = open(args['outFile'], "w")
    if os.path.isfile(args['required_region']):
        log.info("utilizing window file: %s" % args['required_region'])
        window_file = args['required_region']
    else:
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
            chr_start = each_bed[0] + "\t" + each_bed[1] + "\t" + each_bed[2]
            if score is not None:
                outBedGraph.write("%s\t%s\n" % (chr_start, score))
            else:
                outBedGraph.write("%s\tnan\n" % (chr_start))
        outBedGraph.close()
        log.info('done!')

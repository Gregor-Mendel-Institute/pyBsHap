"""
Pakaged functions to get the methyaltion levels using bed tools
"""
import os.path
import numpy as np
from subprocess import Popen, PIPE
import logging
import pandas as pd
import pybedtools as pybed

log = logging.getLogger(__name__)
from pygenome import genome
tair10 = genome.GenomeClass("at_tair10")

def sort_bed_df(bed_df):
    new_bed_df = bed_df.copy()
    new_bed_df = new_bed_df.rename(columns={0: "chr", 1: "start", 2: "end"})
    new_bed_df['strand'] = "+"
    new_bed_df.loc[new_bed_df['start'] > new_bed_df['end'],'strand'] = '-'
    new_bed_df.loc[new_bed_df['strand'] == '-', ['start','end']] = new_bed_df.loc[new_bed_df['strand'] == '-', ['end','start']].values
    new_bed_df = new_bed_df.sort_values(['chr', 'start'])
    return(new_bed_df)

def identify_positions_given_names(in_file, araport11_file):
    if araport11_file is None:
        raise NameError("please give araport file to get positions")
    araport11 = pd.read_csv(araport11_file, header = None, sep = "\t")
    bed_names = np.array(pd.read_csv(in_file, header = None)[0])
    req_bed_df = araport11.loc[araport11[3].isin(bed_names),]
    return(req_bed_df)

def intersect_positions_bed(reference_bed, query_bed):
    assert isinstance(query_bed, pd.DataFrame), "provide a dataframe object" 
    assert query_bed.shape[1] == 2, "provide a dataframe object with only two columns, else use `get_intersect_bed_ix` function"
    assert isinstance(reference_bed, pd.DataFrame), "provide a dataframe object"
    query_ix = np.zeros(0, dtype = "int")
    # ref_ix = np.zeros(0, dtype = "int")
    common_chrs = np.intersect1d(reference_bed.iloc[:,0].unique(), query_bed.iloc[:,0].unique() )
    assert len(common_chrs) > 0, "none of the chromosome IDs are same between reference and query"
    for e_chr in common_chrs:
        e_ref_chr_ix = np.arange( 
            np.searchsorted(reference_bed.iloc[:,0], e_chr, 'left'  ), 
            np.searchsorted(reference_bed.iloc[:,0], e_chr, 'right')
        )
        e_query_chr_ix = np.arange( 
            np.searchsorted(query_bed.iloc[:,0], e_chr, 'left'  ), 
            np.searchsorted(query_bed.iloc[:,0], e_chr, 'right')
        )
        e_ref_chr_pos = np.sort( np.concatenate(reference_bed.iloc[e_ref_chr_ix,:].apply(lambda x: np.arange(x.iloc[1], x.iloc[2] + 1), axis = 1 ).values).ravel() )
        e_query_chr_pos = query_bed.iloc[e_query_chr_ix,1].values

        e_intersect_pos = np.intersect1d(ar1 = e_ref_chr_pos, ar2 = e_query_chr_pos, return_indices = True)
        query_ix = np.append(query_ix, e_query_chr_ix[0] + e_intersect_pos[2])
    return(np.sort(query_ix))

def intersect_positions_bed_minimal(reference_bed_df, query_bed_df, genome_class):
    assert isinstance(query_bed_df, pd.DataFrame), "provide a dataframe object" 
    assert query_bed_df.shape[1] == 2, "provide a dataframe object with only two columns, else use `get_intersect_bed_ix` function"
    assert isinstance(reference_bed_df, pd.DataFrame), "provide a dataframe object"
    assert type(genome_class) is genome.GenomeClass, "Provide a genome class from pygenome repository"

    reference_bed = reference_bed_df.reset_index(drop=True).copy()
    query_bed = query_bed_df.reset_index(drop=True).copy()

    query_genomewide_ix = genome_class.get_genomewide_inds( query_bed )
    ref_genomewide_start_ix = genome_class.get_genomewide_inds( reference_bed.iloc[:,[0,1]] )
    ref_genomewide_ix = np.zeros(0, dtype = "int")
    ref_bed_ix = np.zeros(0, dtype = "int")
    
    reference_bed = sort_bed_df(reference_bed)

    for ef_ref in reference_bed.iterrows():
        ef_ref_pos = np.arange(0, ef_ref[1].iloc[2] - ef_ref[1].iloc[1] + 1 ).astype(int)
        ref_genomewide_ix = np.append(ref_genomewide_ix, ref_genomewide_start_ix[ef_ref[0]] + ef_ref_pos )
        ref_bed_ix = np.append(ref_bed_ix, np.repeat( ef_ref[0], ef_ref_pos.shape[0] )  )
        
    intersect_bed = np.intersect1d(ref_genomewide_ix, query_genomewide_ix, return_indices = True )
    return( (ref_bed_ix[intersect_bed[1]], intersect_bed[2] ) )


def get_intersect_bed_ix(reference_bed, query_bed, just_names=True, araport11_file=None):
    ## here query_bed is either a file or a pandas dataframe
    ## we can rely on bedops -- very fast and efficient
    # https://www.biostars.org/p/319840/
    if isinstance(query_bed, str):
        if os.path.isfile(query_bed):
            queryBed = pybed.BedTool(query_bed)
    elif isinstance(query_bed, pd.DataFrame):
        queryBed = pybed.BedTool.from_dataframe(query_bed.iloc[:,[0,1,2]])
    elif isinstance(query_bed, pybed.bedtool.BedTool):
        queryBed = query_bed
    else:
        raise(NotImplementedError("either input a bed file or pandas dataframe for query"))
    if isinstance(reference_bed, str):
        if os.path.isfile(reference_bed):
            refBed = pybed.BedTool(reference_bed)
    elif just_names:
        reference_bed_df = identify_positions_given_names(reference_bed, araport11_file)
        refBed = pybed.BedTool.from_dataframe(reference_bed_df.iloc[:,[0,1,2]])
    elif isinstance(reference_bed, pd.DataFrame):
        refBed = pybed.BedTool.from_dataframe(reference_bed.iloc[:,[0,1,2]])
    elif isinstance(reference_bed, pybed.bedtool.BedTool):
        refBed = reference_bed
    else:
        raise(NotImplementedError("either input a bed file or pandas dataframe for reference"))
    f_newrefBed = open( refBed.fn + ".new.tmp", 'w' )
    cmd_out = Popen( ''' awk '{ print $0 "\t" NR-1 }' ''' + refBed.fn, shell=True, stdout = f_newrefBed)
    cmd_out.wait()
    f_newrefBed.close()
    newRefBed = pybed.BedTool( refBed.fn + ".new.tmp" )
    f_newqueryBed = open( queryBed.fn + ".new.tmp", 'w' )
    cmd_out = Popen( ''' awk '{ print $0 "\t" NR-1 }' ''' + queryBed.fn, shell=True, stdout = f_newqueryBed)
    cmd_out.wait()
    f_newqueryBed.close()
    newqueryBed = pybed.BedTool( queryBed.fn + ".new.tmp" )
    ## Just taking first three columns for bedtools
    unionBed = newRefBed.intersect(newqueryBed, wa=True, wb = True)
    if unionBed.count() == 0:   ## Return if there are no matching lines.
        return(None)
    unionBed = unionBed.to_dataframe()
    unionBed.columns = np.array(['ref_chr', 'ref_start', 'ref_end', 'ref_ix', 'query_chr', 'query_start', 'query_end', 'query_ix'])
    return(unionBed) ## third column is the index I added 
    # refBed_df = refBed.to_dataframe() 
    # refBed_ids = np.array(refBed_df.iloc[:,0].astype(str) + "," + refBed_df.iloc[:,1].astype(str) + "," +  refBed_df.iloc[:,2].astype(str), dtype="str")
    # unionBed_df = unionBed.to_dataframe() ## wa is to return the entire bed.
    # unionBed_ids = np.array(unionBed_df.iloc[:,0].astype(str) + "," + unionBed_df.iloc[:,1].astype(str) + "," +  unionBed_df.iloc[:,2].astype(str), dtype="str")
    # return(np.where(np.in1d(refBed_ids, unionBed_ids))[0])

## deprecated use the above function for all the chromosomes.
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
        filter_inds.extend(list(range(pos_inds[0],pos_inds[1])))
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
    for echr, echrlen in zip(tair10.chrs, tair10.golden_chrlen):
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
    assert len(feature) < 6
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
    for echr, echrlen in zip(tair10.chrs, tair10.golden_chrlen):
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

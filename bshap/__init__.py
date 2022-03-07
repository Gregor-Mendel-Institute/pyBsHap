"""
    pyBsHap
    ~~~~~~~~~~~~~
    The main module for running pyBsHap
    :copyright: year by my name, see AUTHORS for more details
    :license: license_name, see LICENSE for more details
"""
import os
import os.path
import argparse
import sys
from bshap.core import prebshap
# from bshap.core import bsseq
from bshap.core import bamEdit
from bshap.core import meth5py
from bshap.core import combinemeths
from bshap.core import the1001g
import logging, logging.config

__version__ = '1.1.0'
__updated__ = "17.01.2018"
__date__ = "10.12.2016"

log = logging.getLogger(__name__)
def setLog(logDebug):
    log = logging.getLogger()
    if logDebug:
        numeric_level = getattr(logging, "DEBUG", None)
    else:
        numeric_level = getattr(logging, "ERROR", None)
    log_format = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    lch = logging.StreamHandler()
    lch.setLevel(numeric_level)
    lch.setFormatter(log_format)
    log.setLevel(numeric_level)
    log.addHandler(lch)

def die(msg):
    sys.stderr.write('Error: ' + msg + '\n')
    sys.exit(1)

def get_options(program_license,program_version_message):
    inOptions = argparse.ArgumentParser(description=program_license)
    inOptions.add_argument('-V', '--version', action='version', version=program_version_message)
    subparsers = inOptions.add_subparsers(title='subcommands',description='Choose a command to run',help='Following commands are supported')

    methbam = subparsers.add_parser('getmeth', help="Get methylation on each read from the aligned bam files")
    methbam.add_argument("-i", "--input_bam", dest="inFile", help="aligned BAM file for bs-seq reads")
    methbam.add_argument("-r", "--fasta-file", dest="fastaFile", help="Reference fasta file, TAIR10 genome")
    methbam.add_argument("-s", "--specificRegion", dest="reqRegion", help="region to be checked, Ex. Chr1,1,100 --- an aln file is generated given this", default = '0,0,0')
    methbam.add_argument("-o", "--output", dest="outFile", help="Output file with the methylation across windows")
    methbam.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
    methbam.set_defaults(func=bshap_methbam)

    mhlparser = subparsers.add_parser('getmhl', help="Get methylation haplotype load from the read data from the aligned bam files")
    mhlparser.add_argument("-i", "--input_bam", dest="inFile", help="aligned BAM file for bs-seq reads")
    mhlparser.add_argument("-r", "--fasta-file", dest="fastaFile", help="Reference fasta file, TAIR10 genome")
    mhlparser.add_argument("-d", "--in_hdf5", dest="inhdf5", help="hdf5 file generated using meth5py", default='')
    mhlparser.add_argument("-w", "--window_size", dest="window_size", help="window size", default = 100, type=int)
    mhlparser.add_argument("-x", "--specificRegion", dest="reqRegion", help="region to be checked, Ex. Chr1,1,100 --- an aln file is generated given this", default = '0,0,0')
    mhlparser.add_argument("-o", "--output", dest="outFile", help="Output file", default="STDOUT")
    mhlparser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
    mhlparser.set_defaults(func=bshap_mhlcalc)

    t1001gparser = subparsers.add_parser('generate_h5_1001g', help="generate hdf5 file for input bed files and values.")
    t1001gparser.add_argument("-i", dest="file_paths", nargs='+', help="input bed files containing float variables. You can also provide paths using bash")
    t1001gparser.add_argument("-d", dest="value_column", default=4, type=int, help="column of the files that has values. only if you provide a non-bed")
    t1001gparser.add_argument("-m", dest="matrix", action="store_true",help="mention this option if the file is a matrix. The format is similar to csvsr of R/qtl")
    t1001gparser.add_argument("-o", dest="output_file", help="output h5 file")
    t1001gparser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
    t1001gparser.set_defaults(func=generate_hdf5)

    meth_h5_p = subparsers.add_parser('allc_to_hdf5', help="Generate h5 file from allc")
    meth_h5_p.add_argument("-i", "--input_file", dest="input_file", help="Input file", required=True)
    meth_h5_p.add_argument("-t", "--file_type", dest="file_type", help='Input file type, currently accepted are "bismark" and "methylpy"', required=True)
    meth_h5_p.add_argument("-f", dest="ref_fasta", type=str, help="Path for reference fasta file")
    meth_h5_p.add_argument("--umeth_control", dest="umeth", help='Unmethylated control to calculate conversion rate', default = None)
    meth_h5_p.add_argument("-o", "--output", dest="output_file", help="output file.")
    meth_h5_p.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
    meth_h5_p.set_defaults(func=write_meth_h5file)

    permeth_parser = subparsers.add_parser('methylation_percentage', help="Get methylation percentage on the given bin position.")
    permeth_parser.add_argument("-i", "--input_file", dest="inFile", help="Input methylation HDF5 file generated from allc files", required=True)
    permeth_parser.add_argument("-a", "--allc_path", dest="allc_path", help="Provide this option when you need to create hdf5 file from allc file. 1) Bash path, allc files from this path with sample ID are used.", default="")
    permeth_parser.add_argument("-b", "--required_region", dest="required_region", help="Bed region to calculate the methylation averages. ex. Chr1,1,100", default = '0,0,0')
    permeth_parser.add_argument("-w", "--window_size", dest="window_size", help="window size to get the methylation averages", type = int, default=200)
    permeth_parser.add_argument("-c", "--methylation_averaging_method", dest="category", help="different methylation average methods, 1 -- weighted average, 2 -- methylation calls, 3 -- average methylation per position", default = 1, type = int)
    permeth_parser.add_argument("-o", "--output", dest="outFile", help="output file.")
    permeth_parser.add_argument("--ref_genome", dest="ref_genome", default="at_tair10", type=str, help="Path for reference fasta file")
    permeth_parser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
    permeth_parser.set_defaults(func=bsseq_meth_average)

    cds_meth_parser = subparsers.add_parser('methylation_percentage_cds', help="Get methylation percentage on CDS for given genes")
    cds_meth_parser.add_argument("-i", "--input_file", dest="inFile", help="Input methylation HDF5 file generated from allc files", required=True)
    cds_meth_parser.add_argument("-g", "--gene_bed", dest="gene_bed", help="Bed file for the genes containing gene ID on 4th column", required = True)
    cds_meth_parser.add_argument("-db", "--gffutils_db", dest="gffutils_db", help="gffutils DB, binary file generated from the package.", required = True )
    cds_meth_parser.add_argument("-o", "--output", dest="outFile", help="output file prefix")
    cds_meth_parser.add_argument("--ref_genome", dest="ref_genome", default="at_tair10", type=str, help="Path for reference fasta file")
    cds_meth_parser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
    cds_meth_parser.set_defaults(func=bsseq_cds_meth_average)

    makeCombiAllc = subparsers.add_parser('mergeallc', help = "Merge input allc files from hdf5 format")
    makeCombiAllc.add_argument("-i", dest = "inFiles", nargs='+', help = "Give a list of allc files (hdf5) in here")
    makeCombiAllc.add_argument("-o", dest = "outFile", help = "output hdf5 file combining all the allc")
    makeCombiAllc.add_argument("--min_mc_total", dest = "read_threshold", type = int, help = "filter positions with a minimum read depth. If depth is less than given threshold the permeth is -1", default = 3)
    makeCombiAllc.add_argument("--num_info_lines", dest = "num_info_lines", type = int, help = "filter positions with #lines informative.", default = 20)
    makeCombiAllc.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
    makeCombiAllc.set_defaults(func=mergeallc)

    modifybam = subparsers.add_parser('modifymdtag', help="Modify MD tag on BAM files to get default coloring in jbrowse")
    modifybam.add_argument("-i", "--input_bam", dest="inFile", help="aligned BAM file for bs-seq reads")
    modifybam.add_argument("-r", "--fasta-file", dest="fastaFile", help="Reference fasta file, TAIR10 genome")
    modifybam.add_argument("-o", "--output", dest="outFile", help="Output bam file with modified md tag")
    modifybam.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
    modifybam.set_defaults(func=bshap_modifybam)

    # callmc_parser = subparsers.add_parser('callmc', help="Call mc using methylpy from fastq file")
    # callmc_parser.add_argument("-i", "--input_file", dest="inFile", help="Input fastq file for methylpy")
    # callmc_parser.add_argument("-s", "--sample_id", dest="sample_id", help="unique sample ID for allc Files")
    # callmc_parser.add_argument("-r", "--ref_fol", dest="ref_fol", help="methylpy reference folder for indices and refid", default="/home/GMI/rahul.pisupati/TAiR10_ARABIDOPSIS/03.methylpy.indices/tair10")
    # callmc_parser.add_argument("-f", "--ref_fasta", dest="ref_fasta", help="reference fasta file", default="/home/GMI/rahul.pisupati/TAiR10_ARABIDOPSIS/TAIR10_wholeGenome.fasta")
    # callmc_parser.add_argument("-n", "--nt", dest="nt", help="number of threads", default=2,type=int)
    # callmc_parser.add_argument("-c", "--unMethylatedControl", dest="unMeth", help="unmethylated control", default="ChrC:")
    # callmc_parser.add_argument("-m", "--mem", dest="memory", help="memory for sorting", default="2G")
    # callmc_parser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
    # callmc_parser.set_defaults(func=callmcs_onesample)

    # dmr_parser = subparsers.add_parser('dmrfind', help="Identify DMR using methylpy")
    # dmr_parser.add_argument("-s", "--sample_ids", dest="sample_ids", help="sample ids, comma seperated")
    # dmr_parser.add_argument("-r", "--sample_categories", dest="sample_cat", help="sample categories indicating replicates, comma separated", default="0")
    # dmr_parser.add_argument("-p", "--path", dest="path_to_allc", help="path to allc files")
    # dmr_parser.add_argument("-c", "--context", dest="mc_type", help="methylation context, context separated")
    # dmr_parser.add_argument("-n", "--nt", dest="nt", help="number of threads", default=2,type=int)
    # dmr_parser.add_argument("-o", "--outDMR", dest="outDMR", help="output file for DMR")
    # dmr_parser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
    # dmr_parser.set_defaults(func=dmrfind)

    # lowfreq_parser = subparsers.add_parser('callLowFreq', help="Get lowfreq positions from allc files")
    # lowfreq_parser.add_argument("-s", "--sample_id", dest="sample_id", help="unique sample ID for allc Files")
    # lowfreq_parser.add_argument("-p", "--path", dest="path_to_allc", help="path to allc files")
    # lowfreq_parser.add_argument("-c", "--unMethylatedControl", dest="unMeth", help="unmethylated control", default="ChrC")
    # lowfreq_parser.add_argument("-e", "--pvalue_thres", dest="pvalue_thres", help="threshold for p-value to call low-freq site", default=0.05)
    # lowfreq_parser.add_argument("-o", "--outFile", dest="outFile", help="output h5py file")
    # lowfreq_parser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
    # lowfreq_parser.set_defaults(func=lowfindfind)

    return(inOptions)

def checkARGs(args):
    if not args['inFile']:
        die("input file not specified")
    if not args['fastaFile']:
        die("fasta file not specified")
    if not args['outFile']:
        die("output file not specified")
    if not os.path.isfile(args['fastaFile']):
        die("fasta file does not exist: " + args['fastaFile'])
    if not os.path.isfile(args['inFile']):
        die("input file does not exist: " + args['inFile'])

def bshap_methbam(args):
    checkARGs(args)
    if os.path.isfile(args['reqRegion']):
        prebshap.getMethsRegions(args['inFile'], args['fastaFile'], args['outFile'], args['reqRegion'])
    else:
        prebshap.getMethGenome(args['inFile'], args['fastaFile'], args['outFile'], args['reqRegion'])

def bshap_mhlcalc(args):
    checkARGs(args)
    prebshap.potato_mhl_calc(args)

def mergeallc(args):
    if len(args['inFiles']) <= 1:
        die("provide a list of paths for hdf5 files")
    log.info("reading input files using meth5py")
    meths = combinemeths.CombinedMethsTable( args['inFiles'], file_ids =  None)
    meths_mcs = meths.derive_most_common_positions( args['outFile'], min_mc_total = args['read_threshold'], num_lines_with_info=args['num_info_lines'])
    log.info("done!")

def write_meth_h5file(args):
    m = meth5py.writeHDF5MethTable( ref_fasta_file = args['ref_fasta'] )
    if args['file_type'] == "allc":
        conv_rate = m.load_allc_file( args['input_file'], umeth = args['umeth'] )
    elif args['file_type'] == 'bismark':
        conv_rate = m.load_bismark_coverage( args['input_file'], umeth = args['umeth'] )
    m.write_h5_file( args['output_file'] + ".hdf5" )
    if conv_rate is not None:
        with open(args['output_file'] + '.conv_rate.txt', "w") as out_conv_rate:
            out_conv_rate.write("%s\t%s\n" % (args['input_file'], conv_rate))


def generate_hdf5(args):
    writeh5 = the1001g.WriteHDF51001Table(args['file_paths'], args['output_file'])
    if args['matrix']:
        writeh5.write_h5_matrix()
    else:
        writeh5.write_h5_multiple_files(value_column = args['value_column'])

def bshap_modifybam(args):
    checkARGs(args)
    bamEdit.potatoskin_modify_mdtag_bam(args['inFile'], args['fastaFile'], args['outFile'])

# def callMPsfromVCF(args):
#     if not args['inFile']:
#         die("input file not specified")
#     if not args['outFile']:
#         die("output file not specified")
#     if not os.path.isfile(args['inFile']):
#         die("input file does not exist: " + args['inFile'])
#     bsseq.getMPsfromVCF(args)

# def callmcs_onesample(args):
#     if not args['inFile']:
#         die("input file not specified")
#     if not os.path.isfile(args['inFile']):
#         die("input file does not exist: " + args['inFile'])
#     bsseq.methylpy_callmcs(args)

# def dmrfind(args):
#     bsseq.methylpy_dmrfind(args)

# def lowfindfind(args):
#     bsseq.getLowFreqSites(args)

def bsseq_meth_average(args):
    meth5py.potatoskin_methylation_averages(args)

def bsseq_cds_meth_average(args):
    meth5py.potatoskin_calculate_gbm_exon_only(args)


def main():
    ''' Command line options '''
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = "The main module for pyBsHap"
    program_license = '''%s
    Created by Rahul Pisupati on %s.
    Copyright 2016 Gregor Mendel Institute. All rights reserved.

    Distributed on an "AS IS" basis without warranties
    or conditions of any kind, either express or implied.
    USAGE
    ''' % (program_shortdesc, str(__date__))
    parser = get_options(program_license,program_version_message)
    args = vars(parser.parse_args())
    setLog(args['logDebug'])
    if 'func' not in args:
        parser.print_help()
        return(0)
    try:
        args['func'](args)
        return(0)
    except KeyboardInterrupt:
        return(0)
    except Exception as e:
        logging.exception(e)
        return(2)

if __name__=='__main__':
    sys.exit(main())

import sys
import logging
import argparse
import os


scriptpath = os.path.dirname(sys.argv[0])
sys.path.append(os.path.join(os.path.abspath(scriptpath), "../bshap/core"))

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


import prebshap

inOptions = argparse.ArgumentParser()
subparsers = inOptions.add_subparsers(title='subcommands',description='Choose a command to run',help='Following commands are supported')
methbam = subparsers.add_parser('getmeth', help="pyBsHap on the bam files")
methbam.add_argument("-i", "--input_bam", dest="inFile", help="aligned BAM file for bs-seq reads")
methbam.add_argument("-r", "--fasta-file", dest="fastaFile", help="Reference fasta file, TAIR10 genome")
methbam.add_argument("-o", "--output", dest="outFile", help="Output file with the methylation across windows")
methbam.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")

args = vars(inOptions.parse_args())
setLog(args['logDebug'])



prebshap.getMethWind(args['inFile'], args['fastaFile'], args['outFile'])

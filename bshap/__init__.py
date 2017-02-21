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
import logging, logging.config

__version__ = '0.0.1'
__updated__ = "12.12.2016"
__date__ = "10.12.2016"

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
  methbam = subparsers.add_parser('getmeth', help="pyBsHap on the bam files")
  methbam.add_argument("-i", "--input_bam", dest="inFile", help="aligned BAM file for bs-seq reads")
  methbam.add_argument("-r", "--fasta-file", dest="fastaFile", help="Reference fasta file, TAIR10 genome")
  methbam.add_argument("-s", "--specificRegion", dest="reqRegion", help="region to be checked, Ex. Chr1,1,100 --- an aln file is generated given this", default = '0,0,0')
  methbam.add_argument("-o", "--output", dest="outFile", help="Output file with the methylation across windows")
  methbam.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
  methbam.set_defaults(func=bshap_methbam)
  return inOptions

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
    prebshap.getMethGenome(args['inFile'], args['fastaFile'], args['outFile'], args['reqRegion'])

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
    return 0
  try:
    args['func'](args)
    return 0
  except KeyboardInterrupt:
    return 0
  except Exception as e:
    logging.exception(e)
    return 2

if __name__=='__main__':
  sys.exit(main())

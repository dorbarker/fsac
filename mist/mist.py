import argparse
import json
import subprocess
import io
import sys
import logging
from datetime import datetime
from pathlib import Path

import pandas as pd

logging.basicConfig(filename='.mist.log', level=logging.DEBUG)

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('-t', '--test',
                        type=Path,
                        required=True,
                        help='Test info file')

    parser.add_argument

    parser.add_argument('-a', '--alleles',
                        type=Path,
                        required=True,
                        help='Alleles directory')

    parser.add_argument('-j', '--json-out',
                        type=Path,
                        default=sys.stdout,
                        help='JSON output filename [-]')

    args = parser.parse_args()

    
    if not args.test.exists():
        msg = 'Test file {} does not exist'.format(args.test)
        logging.error(msg)
        raise IOError(msg)

    if not args.allele_dir.exists():
        msg = 'Allele directory {} does not exist'.format(args.alleles)
        logging.error(msg)
        raise IOError(msg)

    return args


def main():

    args = arguments()


def blast(query_path: Path, genome_path: Path) -> pd.DataFrame:

    out_format = '"6 qseqid sseqid pident length qstart qend sstart send qlen \
                  bitscore gaps sseq qseq mismatch"'
    
    blast_cmd = ('blastn',
                 '-query', query_path,
                 '-subject', genome_path,
                 '-outfmt', out_format)

    result = subprocess.run(blast_cmd, check=True, stdout=subprocess.PIPE)

    blast_output = pd.read_csv(io.BytesIO(result.stdout), sep='\t', header=0)

    blast_output['reverse_complement'] = blast_output.sstart > blast_output.send
    
    temp = blast_output[blast_output.reverse_complement]

    blast_output.loc[blast_output.reverse_complement, 'send'] = temp.sstart

    blast_output.loc[blast_output.reverse_complement, 'sstart'] = temp.send

    return blast_output


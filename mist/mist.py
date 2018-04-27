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

    out_format = '6 qseqid sseqid pident length qstart qend sstart send qlen \
                  slen bitscore gaps sseq qseq mismatch'.split()
    
    blast_cmd = ['blastn',
                 '-query', query_path,
                 '-subject', genome_path,
                 '-outfmt'] + out_format

    result = subprocess.run(blast_cmd, check=True, stdout=subprocess.PIPE)

    blast_output = pd.read_csv(io.BytesIO(result.stdout), sep='\t', header=0)

    return blast_output


def parse_blast_results(blast_output: pd.DataFrame) -> pd.DataFrame:

    def is_contig_truncation(row):

        if not row['qlen'] > row['length']:

            return False

        if row['reverse_complement']:

            return row['sstart'] == row['slen'] or row['send'] == 1

        else:

            return row['slen'] == row['send'] or row['sstart'] == 1

    def is_correct(row):

        _correct = (row['mismatch'] == 0) & \
                   (row['pident'] == 100) & \
                   (row['qlen'] == row['slen']) & \
                   (row['gaps'] == 0)

        return _correct

    blast_output['reverse_complement'] = blast_output.sstart > blast_output.send

    correct = blast_output.apply(is_correct, axis=1)

    blast_output['correct'] = correct

    is_truncation = blast_output.apply(is_contig_truncation, axis=1)

    blast_output['is_contig_truncation'] = is_truncation

    return blast_output


def filter_results(blast_output: pd.DataFrame) -> pd.DataFrame:

    # Perfect match
    if any(blast_output['correct']):

        return blast_output[blast_output['correct']]

    # Highest bitscore
    else:

        highest_bitscore = max(blast_output['bitscore'])

        return blast_output[blast_output['bitscore'] == highest_bitscore]


def json_convert(best_blast_hits: pd.DataFrame) -> None:
    pass

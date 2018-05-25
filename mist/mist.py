import argparse
import json
import subprocess
import io
import sys
import re
import logging
from pathlib import Path
from typing import List
import pandas as pd

logging.basicConfig(filename='.mist.log', level=logging.DEBUG)


def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input',
                        type=Path,
                        required=True,
                        help='Input genome')

    parser.add_argument('-t', '--test',
                        type=Path,
                        required=True,
                        help='Test info file')

    parser.add_argument('-a', '--alleles',
                        # nargs='+',
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

    if not args.alleles.exists():
        msg = 'Allele directory {} does not exist'.format(args.alleles)
        logging.error(msg)
        raise IOError(msg)

    try:
        if args.compatible:
            raise NotImplementedError

    except NotImplementedError as e:
        print(e, file=sys.stderr)

    return args


def main():

    args = arguments()

    genome = args.input

    genes = list(args.alleles.glob('*.fasta'))

    blast_results = get_blast_results(genes, genome)

    json_convert(genes, blast_results, args.json_out)


def blast(query_path: Path, genome_path: Path) -> pd.DataFrame:

    out_format = '6 qseqid sseqid pident length qstart qend sstart send qlen \
                  slen bitscore gaps sseq qseq mismatch'
    
    blast_cmd = ['blastn',
                 '-query', str(query_path),
                 '-subject', str(genome_path),
                 '-outfmt', out_format]

    result = subprocess.run(blast_cmd, check=True, stdout=subprocess.PIPE)

    blast_output = pd.read_csv(io.BytesIO(result.stdout),
                               sep='\t',
                               names=out_format.split()[1:])

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
                   (row['qlen'] == row['length']) & \
                   (row['gaps'] == 0)

        return _correct
    # print(max(blast_output['pident']))
    revcomp = blast_output['sstart'] > blast_output['send']
    blast_output['reverse_complement'] = revcomp

    correct = blast_output.apply(is_correct, axis=1)
    blast_output['correct'] = correct

    is_truncation = blast_output.apply(is_contig_truncation, axis=1)
    blast_output['is_contig_truncation'] = is_truncation

    return blast_output


def filter_result(blast_output: pd.DataFrame) -> pd.DataFrame:

    # Perfect match
    if any(blast_output['correct']):

        best = blast_output[blast_output['correct']]

    # Highest bitscore
    else:

        highest_bitscore = max(blast_output['bitscore'])
        best = blast_output[blast_output['bitscore'] == highest_bitscore]

    longest = best[best['length'] == max(best['length'])]

    return longest


def get_blast_results(genes: List[Path], genome: Path) -> List[pd.DataFrame]:

    blast_results = (blast(gene, genome) for gene in genes)

    parsed_results = (parse_blast_results(blast_out)
                      for blast_out in blast_results)

    filtered_results = [filter_result(result) for result in parsed_results]

    return filtered_results


def json_convert(genes: List[Path], best_blast_hits: List[pd.DataFrame],
                 json_out: Path) -> None:

    def marker_match(hit):

        if hit['correct'].item():
            return re.sub('\D*(\d+)', '\\1', hit['qseqid'].item())
        return None

    results = {}

    def unpack_value(value):

        try:
            out = value.item()

        except AttributeError:
            out = value

        return out

    for gene_path, blast_hit in zip(genes, best_blast_hits):

        gene = gene_path.stem

        result = {
            'Mismatches':           blast_hit['mismatch'],
            'QueryAln':             blast_hit['qseq'],
            'SubjAln':              blast_hit['sseq'],
            'Gaps':                 blast_hit['gaps'],
            'QueryName':            blast_hit['qseqid'],
            'SubjName':             blast_hit['sseqid'],
            'PercentIdentity':      blast_hit['pident'],
            'SubjectStartIndex':    blast_hit['sstart'],
            'SubjectEndIndex':      blast_hit['send'],
            'QueryStartIndex':      blast_hit['qstart'],
            'QueryEndIndex':        blast_hit['qend'],
            'BitScore':             blast_hit['bitscore'],
            'ReverseComplement':    blast_hit['reverse_complement'],
            'IsContigTruncation':   blast_hit['is_contig_truncation'],
            'MarkerMatch':          marker_match(blast_hit),
            'CorrectMarkerMatch':   blast_hit['correct'],
        }

        results[gene] = {k: unpack_value(v) for k, v in result.items()}

    with json_out.open('w') as j:
        json.dump(results, j, indent=4)


if __name__ == '__main__':
    main()

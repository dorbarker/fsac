import argparse
import itertools
import json
import logging
import os
import sys
from pathlib import Path

from . import __version__
from .allele_call import allele_call
from .update import update_directory, get_known_alleles
from .tabulate import tabulate_calls

# Ensure numpy, via pandas, doesn't use more than 1 thread.
# If numpy uses multiple threads, it brings no performance benefit in this
# case, but can cause problems if you're running lots of jobs
# on a single HPC node
for env_var in ('OPENBLAS_NUM_THREADS',
                'OMP_NUM_THREADS',
                'MKL_NUM_THREADS',
                'NUMEXPR_NUM_THREADS'):
    os.environ[env_var] = '1'

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('-v', '--version',
                        action='version',
                        version=f'{parser.prog} {__version__}')

    parser.set_defaults(func=None)
    subparsers = parser.add_subparsers(title='Commands')

    ### Allele Calling ###
    allelecall = subparsers.add_parser('call',
                                       help='Call MLST alleles')

    allelecall.set_defaults(func=call_alleles)

    allelecall.add_argument('-i', '--input',
                            type=Path,
                            required=True,
                            help='Input genome')

    allelecall.add_argument('-a', '--alleles',
                            type=Path,
                            required=True,
                            help='Alleles directory')

    allelecall.add_argument('-o', '--output',
                            type=Path,
                            default=sys.stdout,
                            help='JSON output filename [-]')


    ### Allele Updating ###
    update = subparsers.add_parser('update',
                                   help='Update allele definitions')

    update.set_defaults(func=update_results)

    update.add_argument('-a', '--alleles',
                        type=Path,
                        required=True,
                        help='Alleles directory')

    update.add_argument('-j', '--json-dir',
                        type=Path,
                        required=True,
                        help='Directory containing JSON result files')

    update.add_argument('-g', '--genome-dir',
                        type=Path,
                        required=True,
                        help='Directory containing FASTA formatted genomes')

    update.add_argument('-t', '--threshold',
                        type=int,
                        default=10,
                        help='If the hit is this number of basepairs or fewer \
                        shorter than expected, attempt to extend the hit [10]')

    ### Tabulation ###
    tabulate = subparsers.add_parser('tabulate',
                                     help='Create a table from JSON results')

    tabulate.set_defaults(func=tabulate_allele_calls)

    tabulate.add_argument('-j', '--json-dir',
                          type=Path,
                          required=True,
                          help='Directory containing JSON result files')

    tabulate.add_argument('-o', '--output',
                          type=Path,
                          default=sys.stdout,
                          help='Output filename [-]')

    tabulate.add_argument('-d', '--delimiter',
                          type=str,
                          default='\t',
                          help='Delimiter character [TAB]')

    args = parser.parse_args()

    if args.func is None:
        parser.print_help()
        sys.exit(0)

    return args


def main():

    args = arguments()

    logging.basicConfig(datefmt='%Y-%m-%d %H:%M',
                        format='%(asctime)s - %(levelname)s: %(message)s',
                        stream=sys.stderr,
                        level=logging.INFO)

    logging.info('Running %s', args.func.__name__)
    args.func(args)


def validate_fasta(fasta_path: Path):

    if not fasta_path.is_file():
        return (1, f"{fasta_path} does not exist")

    try:
        get_known_alleles(fasta_path)

    except UnboundLocalError:
        return (1, f"{fasta_path} is not in FASTA format")

    except UnicodeDecodeError:
        return (1, f"{fasta_path} is not in FASTA format")

    return (0, "")


def validate_json(json_path: Path):

    if not json_path.is_file():
        return (1, f"{json_path} does not exist")
    try:
        with json_path.open("r") as f:
            data = json.load(f)
            return (0, "")

    except json.decoder.JSONDecodeError:
        return (1, f"{json_path} is not a valid JSON file")


def validate_directory(dir_path: Path, validation_method):

    if not dir_path.is_dir():
        return [(1, f"{dir_path} is not a directory")]

    results = [validation_method(p) for p in dir_path.glob("*")]

    return results


def validate(*args):

    errors, messages = zip(*itertools.chain(*args))

    n_errors = sum(errors)

    if n_errors > 0:

        print(f"Got {n_errors} input errors:")
        print('\n'.join(filter(None, messages)))
        print("Exiting.")
        sys.exit(n_errors)


def call_alleles(args):

    validate([validate_fasta(args.input)],
             validate_directory(args.alleles, validate_fasta))

    allele_call(args.input, args.alleles, args.output)


def update_results(args):

    validate(validate_directory(args.json_dir, validate_json),
             validate_directory(args.alleles, validate_fasta),
             validate_directory(args.genome_dir, validate_fasta))

    update_directory(args.json_dir, args.alleles,
                     args.threshold, args.genome_dir)


def tabulate_allele_calls(args):

    validate(validate_directory(args.json_dir, validate_json))

    tabulate_calls(args.json_dir, args.output, args.delimiter)


if __name__ == '__main__':
    main()

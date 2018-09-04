import argparse
import sys
from pathlib import Path

from src.allele_call import allele_call
from src.update import update_directory
from src.tabulate import tabulate_calls


def arguments():

    parser = argparse.ArgumentParser()

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

    allelecall.add_argument('-j', '--json-out',
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

    args.func(args)


def call_alleles(args):

    allele_call(args.input, args.alleles, args.json_out)


def update_results(args):

    update_directory(args.json_dir, args.alleles)


def tabulate_allele_calls(args):

    tabulate_calls(args.json_dir, args.output, args.delimiter)


if __name__ == '__main__':
    main()

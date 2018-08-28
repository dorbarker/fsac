import argparse
import sys
from pathlib import Path

from src.allele_call import allele_call


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

    allelecall.add_argument('-t', '--test',
                            type=Path,
                            required=True,
                            help='Test info file')

    allelecall.add_argument('-a', '--alleles',
                            type=Path,
                            required=True,
                            help='Alleles directory')

    allelecall.add_argument('-j', '--json-out',
                            type=Path,
                            default=sys.stdout,
                            help='JSON output filename [-]')

    args = parser.parse_args()

    if args.func is None:
        parser.print_help()
        sys.exit(0)

    return args


def main():

    args = arguments()

    args.func(args)


def call_alleles(args):

    allele_call(args.input, args.genes, args.json_out)


if __name__ == '__main__':
    main()

from typing import Dict, Optional, Tuple, Union
from pathlib import Path
import json


def update_locus(gene: Dict[str, Union[str, int, bool, float]],
                 known_alleles: Dict[str, str]) -> Optional[Tuple[str, str]]:
    """
    Modifies known_alleles in place with if a new full-length
    allele has been discovered

    :param gene: Dictionary containing the results of allele_call.allele_call
    :param known_alleles: Dictionary of alleles - sequences are keys,
                                                  headers  are values
    :return: Tuple of sequence and new allele designation
    """

    last_allele = sorted(known_alleles.values(), key=int)[-1]

    next_allele = str(int(last_allele) + 1)

    if gene['CorrectMarkerMatch'] or gene['IsContigTruncation']:

        return None

    seq = gene['SubjAln'].replace('-', '')

    known_alleles[seq] = next_allele

    return seq, next_allele


def update_genome():

    pass


def update_directory():
    pass


def get_known_alleles(alleles_fasta: Path) -> Dict[str, str]:

    known_alleles = {}

    with alleles_fasta.open('r') as f:

        for line in f:

            current_line = line.strip()

            if current_line.strip().startswith('>'):

                try:

                    sequence = ''.join(current_sequence)

                    known_alleles[sequence] = current_header

                except NameError:

                    pass

                current_header = current_line.lstrip('>')

                current_sequence = []

            elif current_line:

                current_sequence.append(current_line)

    return known_alleles

def reverse_complement(sequence: str) -> str:

    complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

    return ''.join([complements[nt] for nt in reversed(sequence)])
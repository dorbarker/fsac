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

    if gene['CorrectMarkerMatch'] or gene['IsContigTruncation']:

        return None, None

    seq = gene['SubjAln'].replace('-', '')

    try:
        allele_name = known_alleles[seq]

    except KeyError:

        last_allele = sorted(known_alleles.values(), key=int)[-1]

        allele_name = str(int(last_allele) + 1)

    return seq, allele_name


def update_genome(genome, genes_dir: Path):

    for gene_name in genome:

        gene = genome[gene_name]

        gene_path = genes_dir / gene_name + '.fasta'

        known_alleles = get_known_alleles(gene_path)

        seq, name = update_locus(gene, known_alleles)

        if seq is None or name is None:
            continue

        # TODO ensure null matches are handled appropriately
        gene['Mismatches'] = 0
        gene['Gaps'] = 0
        gene['QueryName'] = name
        gene['PercentIdentity'] = 100
        gene['MarkerMatch'] = name
        gene['CorrectMarkerMatch'] = True

        genome[gene_name] = gene

def update_directory(results_dir: Path, genes_dir: Path):

    for genome in results_dir.glob('*.json'):

        update_genome(genome, genes_dir)


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
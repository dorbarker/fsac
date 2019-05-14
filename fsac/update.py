from typing import Dict, Optional, Tuple, Union
from pathlib import Path
from contextlib import suppress
import json

SeqAllele = Tuple[Optional[str], Optional[str]]
LocusData = Union[str, int, float, bool]
GeneData = Dict[str, LocusData]


def update_locus(gene: Dict[str, Union[str, int, bool, float]],
                 known_alleles: Dict[str, str]) -> SeqAllele:
    """
    Modifies known_alleles in place with if a new full-length
    allele has been discovered

    :param gene: Dictionary containing the results of allele_call.allele_call
    :param known_alleles: Dictionary of alleles - sequences are keys,
                                                  headers  are values
    :return: Tuple of sequence and new allele designation
    """

    # Gene is missing  - nothing to be done
    if not gene['BlastResult']:
        return None, None

    # Already correct  - nothing to be done
    if gene['CorrectMarkerMatch']:
        return None, None

    # Contig Truncated - nothing to be done
    if gene['IsContigTruncation']:
        return None, None

    # Non-contig trucation
    if gene['PercentLength'] < 1:
        return gene['SubjAln'], None

    seq = gene['SubjAln'].replace('-', '')

    try:
        allele_name = known_alleles[seq]

    except KeyError:

        last_allele = sorted(known_alleles.values(), key=int)[-1]

        allele_name = str(int(last_allele) + 1)

    return seq, allele_name


def update_genome(genome_data: Dict[str, GeneData], genes_dir: Path):

    for gene_name in genome_data:

        gene = genome_data[gene_name]

        gene_path = (genes_dir / gene_name).with_suffix('.fasta')

        known_alleles = get_known_alleles(gene_path)

        seq, name = update_locus(gene, known_alleles)

        if seq is None and name is None:
            continue

        # If there is a non-contig truncation
        if seq is not None and name is None:

            gene['Partial'] = True
            continue

        # TODO ensure null matches are handled appropriately
        gene['Mismatches'] = 0
        gene['Gaps'] = 0
        gene['QueryName'] = name
        gene['PercentIdentity'] = 100
        gene['MarkerMatch'] = name
        gene['CorrectMarkerMatch'] = True

        genome_data[gene_name] = gene

        if seq not in known_alleles:

            update_known_alleles(name, seq, gene_path)

    return genome_data


def update_directory(results_dir: Path, genes_dir: Path):
    """

    :param results_dir: Directory containing JSON results from fsac
    :param genes_dir: Directory containing FASTA files input to fsac
    :return: Void; updates dictionaries in place
    """
    for genome in results_dir.glob('*.json'):

        with genome.open('r') as f:

            genome_data = json.load(f)

            update_genome(genome_data, genes_dir)

        with genome.open('w') as o:

            json.dump(genome_data, o, indent=4)


def get_known_alleles(alleles_fasta: Path) -> Dict[str, str]:

    known_alleles = {}

    with alleles_fasta.open('r') as f:

        for line in f:

            current_line = line.strip()

            if current_line.startswith('>'):

                with suppress(NameError):

                    sequence = ''.join(current_sequence)

                    known_alleles[sequence] = current_header

                current_header = current_line.lstrip('>')

                current_sequence = []

            elif current_line:

                current_sequence.append(current_line)

    return known_alleles


def update_known_alleles(allele_name: str, sequence: str, fasta: Path) -> None:

    with fasta.open('a') as f:

        record = '\n>{name}\n{sequence}'.format(name=allele_name,
                                                sequence=sequence)

        f.write(record)

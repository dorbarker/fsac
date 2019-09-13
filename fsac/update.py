from typing import Dict, Optional, Tuple, Union
from pathlib import Path
from contextlib import suppress
import json

SeqAllele = Tuple[Optional[str], Optional[str]]
LocusData = Union[str, int, float, bool]
GeneData = Dict[str, LocusData]


def update_locus(gene: GeneData,
                 known_alleles: Dict[str, str],
                 threshold: int,
                 genome_path: Path) -> SeqAllele:
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
        seq, _ = extend_hit(gene, threshold, genome_path)

        if seq is None:
            return None, None

    # Full-length, previously unobserved allele
    else:
        seq = gene['SubjAln'].replace('-', '')

    try:
        allele_name = known_alleles[seq]

    except KeyError:

        last_allele = sorted(known_alleles.values(), key=int)[-1]

        allele_name = str(int(last_allele) + 1)

    return seq, allele_name


def extend_hit(gene, threshold: int, genome_path: Path):
    """Extend a BLAST hit if the alignment is less than the threshold shy of
    the expected length. In some cases, it seems that a mismatch near the end
    of the alignment causes the alignment to not be extended.
    """

    difference = gene['QueryLength'] - len(gene['SubjAln'])

    if difference is 0:
        # handle a complete hit
        # return early
        return gene['SubjAln'], gene['MarkerMatch']

    if difference > threshold:
        # handle large discrepancy
        # return early
        return None, None

    # Open subject FASTA
    sequences_names = get_known_alleles(genome_path)
    names_sequences = {value: key for key, value in sequences_names.items()}
    # Find correct contig
    contig = names_sequences[gene['sseqid']]

    # Return target_sequence
    start = gene['sstart'] - 1
    end = gene['sstart'] + gene['qlen']

    full_sequence = contig[start : end]

    if gene['ssend'] > len(contig):
        # handle contig truncation
        return None, None

    assert gene['sseq'] in full_sequence
    return full_sequence, None


def update_genome(genome_data: Dict[str, GeneData],
                  genes_dir: Path,
                  threshold: int,
                  genome_path: Path) -> None:

    for gene_name in genome_data:

        gene = genome_data[gene_name]

        gene_path = (genes_dir / gene_name).with_suffix('.fasta')

        known_alleles = get_known_alleles(gene_path)

        seq, name = update_locus(gene, known_alleles, threshold, genome_path)

        if seq is None and name is None:
            continue

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


def update_directory(results_dir: Path,
                     genes_dir: Path,
                     threshold: int,
                     genomes_path: Path):
    """

    :param results_dir: Directory containing JSON results from fsac
    :param genes_dir: Directory containing FASTA files input to fsac
    :return: Void; updates dictionaries in place
    """
    for genome in results_dir.glob('*.json'):

        genome_path = genomes_path.joinpath(genome.with_suffix('.fasta'))

        with genome.open('r') as f:

            genome_data = json.load(f)

            update_genome(genome_data, genes_dir, threshold, genome_path)

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

        else:

            sequence = ''.join(current_sequence)

            known_alleles[sequence] = current_header

    return known_alleles


def update_known_alleles(allele_name: str, sequence: str, fasta: Path) -> None:

    with fasta.open('a') as f:

        record = '\n>{name}\n{sequence}'.format(name=allele_name,
                                                sequence=sequence)

        f.write(record)

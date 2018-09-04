import json
import pandas as pd
from pathlib import Path
from typing import Dict


def tabulate_calls(jsondir: Path, output: Path, delimiter: str) -> None:

    genomes_calls = compose_table(jsondir)

    write_table(genomes_calls, output, delimiter)


def compose_table(jsondir: Path):

    jsons = jsondir.glob("*.json")

    genomes_calls = {}

    for j in jsons:

        name = j.stem

        calls = dict(parse_gene_calls_from_json(j))

        genomes_calls[name] = calls

    return genomes_calls


def parse_gene_calls_from_json(jsonfile: Path):
    
    with jsonfile.open('r') as f:
        
        for gene in json.load(f):
            
            if not gene['BlastResult']:
                
                call = '0'
                
            elif gene['IsContigTruncation']:
                
                call = '-1'
                
            elif not gene['CorrectMarkerMatch']:
                
                call = '?'
                
            else:
                
                call = gene['MarkerMatch']

            yield gene, call


def write_table(genomes_calls: Dict[str, Dict[str, str]], output: Path,
                delimiter: str) -> None:

    calls = pd.DataFrame(genomes_calls).T

    calls.to_csv(path_or_buf=output, sep=delimiter)

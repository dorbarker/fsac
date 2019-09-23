A Python reimplementation of [Microbial in Silico Typer](https://bitbucket.org/peter87/mist)

To install:
```
conda install -c bioconda -c dorbarker fsac
```

```bash
$ fsac --help
usage: fsac [-h] [-v] {call,update,tabulate} ...

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Print version and exit

Commands:
  {call,update,tabulate}
    call                Call MLST alleles
    update              Update allele definitions
    tabulate            Create a table from JSON results
```

```bash
$ fsac call --help
usage: fsac call [-h] -i INPUT -a ALLELES [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input genome
  -a ALLELES, --alleles ALLELES
                        Alleles directory
  -o OUTPUT, --output OUTPUT
                        JSON output filename [-]
```

```bash
$ fsac update --help
usage: fsac update [-h] -a ALLELES -j JSON_DIR -g GENOME_DIR [-t THRESHOLD]

optional arguments:
  -h, --help            show this help message and exit
  -a ALLELES, --alleles ALLELES
                        Alleles directory
  -j JSON_DIR, --json-dir JSON_DIR
                        Directory containing JSON result files
  -g GENOME_DIR, --genome-dir GENOME_DIR
                        Directory containing FASTA formatted genomes
  -t THRESHOLD, --threshold THRESHOLD
                        If the hit is this number of basepairs or fewer
                        shorter than expected, attempt to extend the hit [10]
```

```bash
$ fsac tabulate --help
usage: fsac tabulate [-h] -j JSON_DIR [-o OUTPUT] [-d DELIMITER]

optional arguments:
  -h, --help            show this help message and exit
  -j JSON_DIR, --json-dir JSON_DIR
                        Directory containing JSON result files
  -o OUTPUT, --output OUTPUT
                        Output filename [-]
  -d DELIMITER, --delimiter DELIMITER
                        Delimiter character [TAB]
```

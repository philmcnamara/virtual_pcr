# Virtual PCR

This program checks off-target binding of primers to the reference genome after design with Primer3 or similar software, with the goal of ensuring the specificity of your PCR

Most in-silico PCR tools only search a predefined set of reference genomes, this tool allows the user the supply any FASTA file as the genome

Output is the position of the expected PCR fragment (chromosome name and base position) along with the expected size

The tool is most easily run by cloning the respository and running the ```virtual_pcr.py``` script from the command line, using the following ```argparse``` flags

## Parameters

| Parameter                | Short Flag | Long Flag     | Note                                                                       |
|--------------------------|------------|---------------|----------------------------------------------------------------------------|
| Reference Genome (FASTA) | -g         | --genome      |                                                                            |
| Upper Size Limit         | -u         | --upper_limit | Maximum size for amplified fragment (bp)                                   |
| Lower Size Limit         | -l         | --lower_limit | Minimum size for amplified fragment (bp)                                   |
| Mismatches               | -m         | --mismatch    | Maximum number of mismatches allowed between primer and genome (default 0) |
| Forward Primer           | -f         | --f_primer    | Forward primer (5' to 3')                                                  |
| Reverse Primer           | -r         | --r_primer    | Reverse primer (5' to 3')                                                  |

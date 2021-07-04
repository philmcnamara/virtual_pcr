from virtual_pcr.virtual_pcr import check_chromosome, generate_chromosome_dict, reverse_complement
from pathlib import Path

def test_reverse_complement():
    assert reverse_complement("ACAGCG") == "CGCTGT"
    assert reverse_complement("GATGATC") == "GATCATC"
    assert reverse_complement("CCTTGGAA") == "TTCCAAGG"

def test_generate_chromosome_dict():
    chrom_dict = generate_chromosome_dict(Path.cwd() / "tests" / "test_genome.fa")
    assert len(chrom_dict) == 2
    # make sure we are reading both single-line and multi-line FASTA sequences correctly
    assert len(chrom_dict["chrom_1"]) == 48
    assert len(chrom_dict["chrom_2"]) == 48

def test_pattern_1_exact_match():
    chrom_dict = generate_chromosome_dict(Path.cwd() / "tests" / "test_genome.fa")
    fake_args = {"f_primer": "GACTGATA", "r_primer": "ATCTGTGT", "mismatch": 0,  "lower_limit": 1, "upper_limit": 100}

    results = check_chromosome(chrom_dict["chrom_1"], fake_args)

    assert len(results["pattern 1"]) == 1 and len(results["pattern 2"]) == 0

def test_mismatch():
    chrom_dict = generate_chromosome_dict(Path.cwd() / "tests" / "test_genome.fa")
    # mismatch - change last nucleotide of f_primer to T
    fake_args = {"f_primer": "GACTGATT", "r_primer": "ATCTGTGT", "mismatch": 0,  "lower_limit": 1, "upper_limit": 100}
    results = check_chromosome(chrom_dict["chrom_1"], fake_args)

    assert len(results["pattern 1"]) == 0 and len(results["pattern 2"]) == 0

    fake_args["mismatch"] = 1
    
    results = check_chromosome(chrom_dict["chrom_1"], fake_args)
    # with 1 allowed mismatch we should get 2 hits, length 25 and 46
    assert len(results["pattern 1"]) == 2 and len(results["pattern 2"]) == 0

def test_pattern_2_match():
    chrom_dict = generate_chromosome_dict(Path.cwd() / "tests" / "test_genome.fa")
    fake_args = {"f_primer": "ATCTGTGTG", "r_primer": "GATTTATAT", "mismatch": 0,  "lower_limit": 1, "upper_limit": 100}
    results = check_chromosome(chrom_dict["chrom_2"], fake_args)
    assert len(results["pattern 2"]) == 1 and len(results["pattern 1"]) == 0
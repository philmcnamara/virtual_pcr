import regex
import argparse

def get_arguments():
    parser = argparse.ArgumentParser(
        description="Virtual PCR Tool"
    )
    parser.add_argument(
        "-g", "--genome",
        help="Reference genome (FASTA)",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-u", "--upper_limit",
        help="Maximum size for amplified fragment (bp)",
        required=True,
        type=int,
    )
    parser.add_argument(
        "-l", "--lower_limit",
        help="Minimum size for amplified fragment (bp)",
        required=True,
        type=int,
    )
    parser.add_argument(
        "-m", "--mismatch",
        help="Maximum number of mismatches allowed between primer and genome (default 3)",
        required=False,
        type=int,
        default=0
    )
    parser.add_argument(
        "-f", "--f_primer",
        help="Forward primer (5' to 3')",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-r", "--r_primer",
        help="Reverse primer (5' to 3')",
        required=True,
        type=str,
    )

    return vars(parser.parse_args())

def reverse_complement(seq):
    base_dict = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join([base_dict[i] for i in seq[::-1]])

def check_chromosome(chrom_seq, in_args):

    results = {}

    # Pattern 1 is the forward primer binding on the forward strand
    # and the reverse primer binding downstream on the reverse strand
    f_primer_1 = in_args["f_primer"]
    r_primer_1 = reverse_complement(in_args["r_primer"])
    results["pattern 1"] = find_matches(chrom_seq, f_primer_1, r_primer_1, in_args)

    # Pattern 2 is the reverse primer binding on the forward strand
    # and the forward primer binding downstream on the reverse strand
    f_primer_2 = in_args["r_primer"]
    r_primer_2 = reverse_complement(in_args["f_primer"])
    results["pattern 2"] = find_matches(chrom_seq, f_primer_2, r_primer_2, in_args)

    return results

def find_matches(chrom_seq, f_primer, r_primer, in_args):

    proper_matches = []

    f_primer_query = "(?:%s){s<=%d}" % (f_primer, in_args["mismatch"])
    r_primer_query = "(?:%s){s<=%d}" % (r_primer, in_args["mismatch"])

    f_matches = regex.finditer(f_primer_query, chrom_seq)
    r_matches = regex.finditer(r_primer_query, chrom_seq)

    # product size is from f start to r stop
    f_starts = [i.start() for i in f_matches]
    r_stops = [i.end() for i in r_matches]

    # check all combinations of f starts and r stops
    for i in f_starts:
        for j in r_stops:
            size = j - i
            if in_args["lower_limit"] <= size <= in_args["upper_limit"]:
                proper_matches.append((i, size))

    return proper_matches

def generate_chromosome_dict(reference_genome):
    chrom_dict = {}
    with open(reference_genome, "r") as genome_in:
        line = genome_in.readline().strip("\n").strip("\r")
        while line:
            if line.startswith(">"):
                chrom_name = line.lstrip(">")
                chrom_list = []
                chrom_line = genome_in.readline().strip("\n").strip("\r")
                while chrom_line and not chrom_line.startswith(">"):
                    chrom_list.append(chrom_line)
                    chrom_line = genome_in.readline().strip("\n").strip("\r")
                line = chrom_line
                chromosome = "".join(chrom_list)
                chrom_dict[chrom_name] = chromosome
    return chrom_dict

def main():
    in_args = get_arguments()

    chrom_dict = generate_chromosome_dict(in_args["genome"])

    for chrom_name in chrom_dict:
        results = check_chromosome(chrom_dict[chrom_name], in_args)
        if len(results["pattern 1"]) > 0 or len(results["pattern 2"]) > 0:
            print("Matches on {}".format(chrom_name))
            print("".join("-" * 50))
            if len(results["pattern 1"]) > 0:
                print("Pattern 1 matches (Forward primer on forward strand)")
                print("".join("-" * 50))
                for hit in results["pattern 1"]:
                    print("Start position: {}".format(hit[0]))
                    print("Length (bp) : {}".format(hit[1]))
                    print("".join("-" * 50))
            if len(results["pattern 2"]) > 0:
                print("Pattern 2 matches (Forward primer on reverse strand")
                for hit in results["pattern 2"]:
                    print("Start position: {}".format(hit[0]))
                    print("Length (bp) : {}".format(hit[1]))
                    print("".join("-" * 50))

if __name__ == "__main__":
   main()

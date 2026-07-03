#!/usr/bin/env python

from __future__ import print_function
import argparse
import gzip
import os


def _read_fasta(fastafile):
    """Read sequences from a FASTA file into a dict keyed by '>header'."""
    sequences = {}
    if fastafile.endswith("gz"):
        with gzip.open(fastafile, 'r') as f:
            for line in f:
                line = str(line, encoding="utf-8")
                if line.startswith(">"):
                    if " " in line:
                        seq, _ = line.split(' ', 1)
                        sequences[seq] = ""
                    else:
                        seq = line.rstrip("\n")
                        sequences[seq] = ""
                else:
                    sequences[seq] += line.rstrip("\n")
    else:
        with open(fastafile, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    if " " in line:
                        seq, _ = line.split(' ', 1)
                        sequences[seq] = ""
                    else:
                        seq = line.rstrip("\n")
                        sequences[seq] = ""
                else:
                    sequences[seq] += line.rstrip("\n")
    return sequences


def _read_cluster_map(resultfile):
    """Read a two-column TSV (contig_name, cluster_name) into a dict."""
    dic = {}
    with open(resultfile, "r") as f:
        for line in f:
            contig_name, cluster_name = line.strip().split('\t')
            dic.setdefault(cluster_name, []).append(contig_name)
    return dic


def gen_bins(fastafile, resultfile, outputdir):
    """Generate bin FASTA files using sequential integer filenames."""
    print("Processing file:\t{}".format(fastafile))
    sequences = _read_fasta(fastafile)
    print("Reading Map:\t{}".format(resultfile))
    dic = _read_cluster_map(resultfile)
    print("Writing bins:\t{}".format(outputdir))
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    bin_name = 0
    for _, cluster in dic.items():
        binfile = os.path.join(outputdir, "{}.fa".format(bin_name))
        with open(binfile, "w") as f:
            for contig_name in cluster:
                contig_name = ">" + contig_name
                try:
                    sequence = sequences[contig_name]
                except KeyError:
                    bin_name += 1
                    continue
                f.write(contig_name + "\n")
                f.write(sequence + "\n")
                bin_name += 1


def gen_bins_with_cluster_ids(fastafile, resultfile, outputdir):
    """Generate bin FASTA files using the cluster ID from the TSV as the filename.

    This preserves the mapping between cluster IDs and bin files, which is
    required to match CheckM2 quality_report.tsv entries back to cluster IDs.

    Parameters
    ----------
    fastafile : str
        Path to the input FASTA file containing contig sequences.
    resultfile : str
        Path to a two-column TSV mapping contig names to cluster IDs.
    outputdir : str
        Directory where per-cluster FASTA files will be written.
    """
    print("Processing file:\t{}".format(fastafile))
    sequences = _read_fasta(fastafile)
    print("Reading Map:\t{}".format(resultfile))
    dic = _read_cluster_map(resultfile)
    print("Writing bins:\t{}".format(outputdir))
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    for cluster_name, cluster in dic.items():
        # Sanitize cluster_name so it is safe to use as a filename
        safe_name = str(cluster_name).replace('/', '_').replace('\\', '_')
        binfile = os.path.join(outputdir, "{}.fa".format(safe_name))
        with open(binfile, "w") as f:
            for contig_name in cluster:
                contig_key = ">" + contig_name
                if contig_key not in sequences:
                    continue
                f.write(contig_key + "\n")
                f.write(sequences[contig_key] + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="original fasta file")
    parser.add_argument("-r", help="tsv version result file")
    parser.add_argument("-o", help="output dir")
    args = parser.parse_args()
    gen_bins(args.f, args.r, args.o)

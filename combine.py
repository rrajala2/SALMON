#!/usr/bin/env python

HELP = """
This program combines kallisto and salmon abundance files scattered across directories.

The output of the program is similar to the count file produced by featureCounts

Usage:
    
    cat IDS | python combine-abundance.py RESULTS 
    
Where:

 - RESULTS is a directory
 - IDS a row oriented file with one identifier per line
 
For each id the IDS file the program will attempt to load the file named as 

   RESULTS/id/abundances.tsv

Once all files are loaded these are combined into a result that has six columns to start each row, 
then the estimated abundance from each file labeled with a sample name.

"""

import csv
import os
import sys


def parse(workdir, sample):
    """
    Returns a generator over a file
    """

    # The path to the abundnace file.
    path = os.path.join(workdir, sample, "abundance.tsv")

    # Try the salmon abundance file.
    if not os.path.isfile(path):
        path = os.path.join(workdir, sample, "quant.sf")

    # Stop on missing files.
    if not os.path.isfile(path):
        print("ERROR: unable to read either kallisto or salmon file: %s" % path)
        sys.exit()

    # Open the stream.
    stream = csv.DictReader(open(path), delimiter="\t")

    # Generate the data
    for row in stream:
        yield row


def process(workdir, stream=sys.stdin):
    """
    Processes and combines tsv files.
    """

    # The samples in the work directory.
    samples = [line.strip() for line in stream]

    collect = {}
    for sample in samples:
        rows = list(parse(workdir, sample))

        for row in rows:

            try:
                target_id, length, eff_length = row['target_id'], row['length'], row['eff_length']
            except KeyError as exc:
                target_id, length, eff_length = row['Name'], row['Length'], row['EffectiveLength']

            # Fill in the first six columns.
            if target_id not in collect:
                collect[target_id] = [target_id, length, eff_length ]

            try:
                value = row['est_counts']
            except KeyError as exc:
                value = row['NumReads']

            value = str(round(float(value), 1))
            collect[target_id].append(value)

    # Print the collected data
    header = [ "target_id", "length", "eff_length" ] + samples
    header = "\t".join(header)
    print (header)

    # Print the header

    for key, values in collect.items():
        line = "\t".join(values)
        print(line)

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print(HELP)
        sys.exit()

    # Get the work directory
    workdir = sys.argv[1]

    if not os.path.isdir(workdir):
        print("ERROR: %s is not a directory" % workdir)
        sys.exit(-1)

    process(workdir=workdir, stream=sys.stdin)

#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "IMBERT Pierre"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["IMBERT Pierre"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "IMBERT Pierre"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file',
                        type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen',
                        type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size',
                        type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


#==============================================================
# Dé-duplication en séquence "complète"
#==============================================================
def read_fasta(amplicon_file, minseqlen):
    """
    Read a .fasta.gz file.

    Parameter
    ---------
    amplicon_file: str
    minseqlen: int

    Return
    ------
    a generator of sequences of size >=  minseqlen
    """
    with gzip.open(amplicon_file, "rb") as filin:
        for _ in filin:
            try:
                seq = next(filin).decode('UTF-8').strip()
                if len(seq) >= minseqlen:
                    yield seq
            except StopIteration:
                return


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """
    Determine unique sequences that occurs at least n time.

    With n >= mincount.

    Parameter
    ---------
    amplicon_file: str
         name of the fasta file
    minseqlen: int
         minimum length of each sequences
    mincount: int
         minimum sequence count

    Return
    ------
    generator of occurrence [sequence, count] - descending order
    """
    occ = {}
    for seq in read_fasta(amplicon_file, minseqlen):
        if not seq in occ:
            occ[seq] = 0
        occ[seq] += 1

    # Sort occ dictionary by value - descending order
    new_occ = {
        k: v for k, v in sorted(occ.items(), key=lambda item: item[1], reverse=True)
    }

    for seq, count in new_occ.items():
        if count >= mincount:
            try:
                yield [seq, count]
            except StopIteration:
                return


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    derep = dereplication_fulllength(args.amplicon_file, args.minseqlen,
                                     args.mincount)


if __name__ == '__main__':
    main()

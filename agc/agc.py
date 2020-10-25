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
    if amplicon_file.endswith(".gz"):
        # Lecture d'un fichier fasta.gz
        with gzip.open(amplicon_file, "rb") as filin:
            for line in filin:
                try:
                    seq = next(filin).decode('UTF-8').strip()
                    if len(seq) >= minseqlen:
                        yield seq
                except StopIteration:
                    return
    else:
        # Lecture d'un fichier fasta
        with open(amplicon_file, "r") as filin:
            seq = ""
            for line in filin:
                if not line.startswith('>'):
                    seq += line.strip()
                else:
                    try:
                        if len(seq) >= minseqlen: yield seq
                    except StopIteration:
                        return
                    seq = ""
            try:
                if len(seq) >= minseqlen: yield seq
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
# Recherche d'une séquence chimérique par approche "de novo"
#==============================================================
def get_chunks(sequence, chunk_size):
    """
    Determine sub-sequences of size I non-overlapping.

    Parameter
    ---------
    sequence: str
    chunk_size: int

    Return
    ------
    segments: list
        list of sub-sequences, at least of size 4
    """
    segments = []
    for i in range(0, len(sequence), chunk_size):
        tmp = sequence[i:chunk_size+i]
        if len(tmp) == chunk_size:
            segments.append(tmp)
    return segments


def cut_kmer(sequence, kmer_size):
    """
    Generate kmer of size k of a sequence.

    Parameter
    ---------
    sequence: str
    kmer_size: int

    Return
    ------
    generator of kmer of size k
    """
    for i in range(len(sequence)-kmer_size+1):
        try:
            yield sequence[i:kmer_size+i]
        except StopIteration:
            return


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """
    Parameter
    ---------
    kmer_dict: dictionary
    sequence: str
    id_seq: int
    kmer_size: int

    Return
    ------
    kmer_dict - key: kmer, value: set of sequence id
    """
    kmer_seq = cut_kmer(sequence, kmer_size)

    for kmer in kmer_seq:
        if not kmer in kmer_dict:
            kmer_dict[kmer] = set()
        kmer_dict[kmer].add(id_seq)

    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
    """
    Parameter
    ---------
    kmer_dict: dictionary
    sequence: str
    kmer_size: int

    Return
    ------
    best_mates: list
    """
    kmer_seq =  cut_kmer(sequence, kmer_size)
    tmp = [
        i[0] for i in Counter([ids for kmer in kmer_seq if kmer in kmer_dict for ids in kmer_dict[kmer]]).most_common(8)
    ]
    return tmp


def get_identity(alignment_list):
    """
    Compute the percent of identity between two sequences.

    identity = n / N avec

      - n: identical nucleotides
      - N: size of alignment

    Parameter
    ---------
    alignment_list: list
    """
    size_alignement = len(alignment_list[1])
    identical_nucleotide = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            identical_nucleotide += 1

    return (identical_nucleotide / size_alignement) * 100


def detect_chimera(perc_identity_matrix):
    """

    Parameter
    ---------
    perc_identity_matrix:
        for each chunk give the identity between a target sequence, that could be a
        chimera, and two parents sequences

    Return
    ------
    True: the sequence is a chimera
    False: the sequence is not a chimera
    """
    chimera = False
    std_chunk = 0
    idt_chunk1, idt_chunk2 = set(), set()  # nombre de segments ayant une identité !=

    for identity in perc_identity_matrix:
        std_chunk += std(identity)
        idt_chunk1.add(identity[0])
        idt_chunk2.add(identity[1])

    if (std_chunk/len(perc_identity_matrix)) > 5 and \
       (len(idt_chunk1) >=2 or len(idt_chunk2) >= 2):
        return True
    return False


def std(values):
    """values: list"""
    return statistics.stdev(values)


def get_unique(ids):
    """Unique."""
    return {}.fromkeys(ids).keys()


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    Determine which sequences are chimera or not.

    Parameter
    ---------
    amplicon_file: str
    minseqlen: int
    mincount: int
    chunk_size: int
    kmer_size: int

    Return
    ------
    generator of non-chimera sequences
    """
    # Sequence
    sequences = []
    occ = []
    for de_rep in dereplication_fulllength(amplicon_file, minseqlen, mincount):
        sequences.append(de_rep[0])
        occ.append(de_rep[1])

    # Séparation en segment de taille chunk_size + génération du dictionnaire de kmer
    segments, kmer_dico = [], {}
    for i in range(len(sequences)):
        segments.append(get_chunks(sequences[i], chunk_size))
        kmer_dico = get_unique_kmer(kmer_dico, sequences[i], i, kmer_size)

    # Génération des best_mates pour un segment donné
    best_mates = []
    for sequence_chunks in segments:
        for each_chunk in sequence_chunks:
            best_mates.append(search_mates(kmer_dico, each_chunk, kmer_size))

    # Recherche de séquences parentes - séquences présentes dans toutes les listes
    seq_parentes = common(best_mates[0], best_mates[1])

    # Déterminer si une séquence est une chimère
    chimera_id = []
    chunk_seq_list = [get_chunks(sequences[seq_parentes[0]], chunk_size)]
    chunk_seq_list += [get_chunks(sequences[seq_parentes[1]], chunk_size)]
    for i in range(len(sequences)):
        if not i in seq_parentes:
            chunk_chim = get_chunks(sequences[i], chunk_size)

            perc_identity_matrix = [[] for c in range(len(chunk_chim))]
            for j in range(len(chunk_seq_list)):
                for l,chunk in enumerate(chunk_chim):
                    perc_identity_matrix[l].append(
                        get_identity(nw.global_align(chunk, chunk_seq_list[j][l], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__), '../agc')) + "/MATCH")))

            if detect_chimera(perc_identity_matrix):
                chimera_id.append(i)


    for i in range(len(sequences)):
        if not i in chimera_id:
            yield [sequences[i], occ[i]]


def common(liste1, liste2):
    """Common."""
    return list(set(liste1) & set(liste2))


#==============================================================
# Regroupement glouton
#==============================================================
def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size,
                                kmer_size):
    """
    Compute the occurrence of each non-chimera sequences.

    Parameter
    ---------
    amplicon_file: str
    minseqlen: int
    mincount: int
    chunk_size: int
    kmer_size: int

    Return
    ------
    otu: list
        count of each sequences
    """
    data = chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)
    otu = []

    for sequence, count in data:
        otu.append((sequence, count))

    return otu


def write_OTU(otu_list, output_file):
    """
    Parameter
    ---------
    otu_list: list
    output_file: str

    Format OTU:
    >OTU_{numéro partant de 1} occurrence:{nombre d’occurrence à la déréplication}
    {séquence au format fasta}
    """
    with open(output_file, "w") as filout:
        for i, otu in enumerate(otu_list):
            filout.write(">OTU_{} occurrence:{}\n".format(i+1, otu[1]))
            filout.write("{}\n".format(fill(otu[0])))


def fill(text, width=80):
    """Split text with a line return to respect fasta format."""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    otu = abundance_greedy_clustering(args.amplicon_file, args.minseqlen,
                                      args.mincount, args.chunk_size, args.kmer_size)

    write_OTU(otu, args.output_file)


if __name__ == '__main__':
    main()

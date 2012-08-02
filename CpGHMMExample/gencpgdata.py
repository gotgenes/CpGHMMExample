#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# Copyright (c) 2012 Christopher D. Lasher
#
# This software is released under the MIT License. Please see
# LICENSE.txt for details.


"""Generates a random sequence with CpG islands.

This script produces three outfiles:

* a FASTA format sequence file

* a file containing the start and end positions of CpG islands

* a file containing the parameters of the transitions
"""

import argparse
import bisect
import random
import textwrap

import numpy as np

import logging
LOGGER = logging.getLogger()
LOGGER.setLevel(logging.INFO)
STREAM_HANDLER = logging.StreamHandler()
STREAM_HANDLER.setLevel(logging.INFO)
LOGGER.addHandler(STREAM_HANDLER)
FORMATTER = logging.Formatter('%(message)s')
STREAM_HANDLER.setFormatter(FORMATTER)

ALPHABET = 'ACGT'

# Rows are ACGT, columns are ACGT
_CPG_CPG_PROBABILITIES = np.array([
    [
        0.1,
        0.4,
        0.4,
        0.1
    ],
    [
        0.05,
        0.45,
        0.45,
        0.05
    ],
    [
        0.05,
        0.45,
        0.45,
        0.05
    ],
    [
        0.1,
        0.4,
        0.4,
        0.1
    ],
])

# Rows are ACGT, columns are ACGT
_NORMAL_NORMAL_PROBABILITIES = np.array([
    [
        0.25,
        0.25,
        0.25,
        0.25
    ],
    [
        0.15,
        0.35,
        0.35,
        0.15
    ],
    [
        0.15,
        0.35,
        0.35,
        0.15
    ],
    [
        0.25,
        0.25,
        0.25,
        0.25
    ],
])

_CPG_TO_NORMAL_TRANSITION_PROB = 0.005
_NORMAL_TO_CPG_TRANSITION_PROB = 0.0025

_CPG_PROBABILITIES = np.concatenate(
    (
        (1 - _CPG_TO_NORMAL_TRANSITION_PROB) * _CPG_CPG_PROBABILITIES,
        _CPG_TO_NORMAL_TRANSITION_PROB * _NORMAL_NORMAL_PROBABILITIES
    ),
    1
)
_NORMAL_PROBABILITIES = np.concatenate(
    (
        _NORMAL_TO_CPG_TRANSITION_PROB * _CPG_CPG_PROBABILITIES,
        (1 - _NORMAL_TO_CPG_TRANSITION_PROB) * _NORMAL_NORMAL_PROBABILITIES
    ),
    1
)

TRANSITION_PROBABILITIES = np.concatenate(
        (_CPG_PROBABILITIES, _NORMAL_PROBABILITIES))

TRANSITION_CUMSUMS = TRANSITION_PROBABILITIES.cumsum(1).tolist()
for row in TRANSITION_CUMSUMS:
    row[-1] = 1.0


def generate_sequence(length):
    """Generates the random sequence, including CpG islands.

    :param length: length of the sequence to generate
    :returns: a randomly generated sequence, and a list of start and end
        positions of CpG sites within the sequence

    """
    sequence = []
    cpg_sites = []
    cpg_start = None
    in_cpg = False
    start = random.randrange(len(TRANSITION_CUMSUMS))
    sequence.append(ALPHABET[start % 4])
    if start < 4:
        in_cpg = True
        cpg_start = start
    prev_index = start
    for x in range(1, length):
        random_value = random.random()
        transition_index = bisect.bisect_left(
                TRANSITION_CUMSUMS[prev_index], random_value)
        sequence.append(ALPHABET[transition_index % 4])
        if transition_index < 4:
            if not in_cpg:
                cpg_start = x
                in_cpg = True
        else:
            if in_cpg:
                cpg_sites.append((cpg_start, x - 1))
                in_cpg = False
        prev_index = transition_index

    if in_cpg:
        cpg_sites.append((cpg_start, length - 1))

    return ''.join(sequence), cpg_sites


def wrap_sequence(sequence, width=50):
    return '\n'.join(sequence[i:i+width] for i in
                     xrange(0, len(sequence), width))


def output_sequence(outfileh, sequence):
    """Writes the sequence to the outfile in FASTA format.

    :param outfileh: @todo
    :param sequence: @todo
    :returns: @todo

    """
    outfileh.write('>testcpg\n')
    formatted_sequence = wrap_sequence(sequence, 50)
    outfileh.write(formatted_sequence)


def output_sites(outfileh, sites):
    """Writes the CpG start and end positions to a CSV-format file.

    :param outfileh: @todo
    :param sites: @todo
    :returns: @todo

    """
    outlines = ("{},{}\n".format(start, end) for (start, end) in sites)
    outfileh.writelines(outlines)


def make_cli_parser():
    """Creates the command-line interface.

    :returns: an :py:class:`argparse.ArgumentParser` instance

    """
    cli_parser = argparse.ArgumentParser(description=__doc__)
    cli_parser.add_argument(
            'length', type=int, help="length of sequence to generate")
    return cli_parser


def main(argv=None):
    cli_parser = make_cli_parser()
    args = cli_parser.parse_args(argv)
    LOGGER.info("Generating random CpG sequence.")
    sequence, cpg_sites = generate_sequence(args.length)
    LOGGER.info("Writing sequence to test_cpg_sequence.fasta")
    with open('test_cpg_sequence.fasta', 'w') as fasta_outfile:
        output_sequence(fasta_outfile, sequence)
    LOGGER.info("Writing CpG site positions to test_cpg_sites.csv")
    with open('test_cpg_sites.csv', 'w') as positions_outfile:
        output_sites(positions_outfile, cpg_sites)
    LOGGER.info("Writing transition probabilities to "
                "test_cpg_transitions.csv")
    np.savetxt('test_cpg_transitions.csv', TRANSITION_PROBABILITIES,
               delimiter=',')


if __name__ == '__main__':
    main()


#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# Copyright (c) 2012 Christopher D. Lasher
#
# This software is released under the MIT License. Please see
# LICENSE.txt for details.


"""Identifies potential CpG sites within a given sequence."""

import argparse
from math import exp, log

import numpy as np
import matplotlib.pyplot as plt

from CpGHMMExample import alglib, gencpgdata

import logging
LOGGER = logging.getLogger()
LOGGER.setLevel(logging.INFO)
#STREAM_HANDLER = logging.StreamHandler()
#STREAM_HANDLER.setLevel(logging.INFO)
#LOGGER.addHandler(STREAM_HANDLER)
#FORMATTER = logging.Formatter('%(message)s')
#STREAM_HANDLER.setFormatter(FORMATTER)

NUC_TO_INDICES = {
    'A': (0, 4),
    'C': (1, 5),
    'G': (2, 6),
    'T': (3, 7)
}


def _parse_fasta_file(fileh):
    fileh.next()
    sequence = ''.join([line.strip() for line in fileh])
    return sequence


def _parse_known_sites_file(fileh):
    sites = []
    for line in fileh:
        start, end = line.split(',')
        sites.append((int(start), int(end)))
    return sites


def log_add(log_value_x, log_value_y):
    """Add two values together in log scale.

    The idea here is that the equation :math:`z = x + y` can be
    re-written as

    .. math::

       \log z = \log x + \log(1 + e^{(\log y - \log y)})

    :param log_value_x: some value x in log space
    :param log_value_y: some value y in log space
    :returns: the sum of ``x`` and ``y`` in log space

    """
    return log_value_x + log(1 + exp(log_value_y - log_value_x))


def calc_forward_probabilities(sequence, transition_matrix):
    """Calculate the forward probabilities for a given sequence and
    transition matrix.

    :param sequence: @todo
    :param transition_matrix: @todo
    :returns: @todo

    """
    num_symbols = transition_matrix.shape[0]
    seq_length = len(sequence)
    forward_arr = np.zeros((num_symbols, seq_length))
    # Starting probability will be considered equal for all states.
    start_prob = log(1 / 8.0)
    prev_symbols = NUC_TO_INDICES[sequence[0]]
    for symbol in prev_symbols:
        forward_arr[symbol, 0] = start_prob
    # Now calculate for remainder of the sequence.
    for i in range(1, seq_length):
        symbols = NUC_TO_INDICES[sequence[i]]
        for symbol in symbols:
            prob_transition_symbol1 = (
                forward_arr[prev_symbols[0], i-1] +
                transition_matrix[symbol, prev_symbols[0]]
            )
            prob_transition_symbol2 = (
                forward_arr[prev_symbols[1], i-1] +
                transition_matrix[symbol, prev_symbols[1]]
            )
            forward_arr[symbol, i] = (
                prob_transition_symbol1 +
                log(1 + exp(prob_transition_symbol2 -
                            prob_transition_symbol1))
            )
        prev_symbols = symbols

    final_transition_prob1 = (
        forward_arr[prev_symbols[0], -1] + start_prob)
    final_transition_prob2 = (
        forward_arr[prev_symbols[1], -1] + start_prob)
    final_probability = (
        final_transition_prob1 +
        log(1 + exp(final_transition_prob2 - final_transition_prob1))
    )
    return forward_arr, final_probability


def calc_backward_probabilities(sequence, transition_matrix):
    """Calculate the backward probabilities for a given sequence and
    transition matrix.

    :param sequence: @todo
    :param transition_matrix: @todo
    :returns: @todo

    """
    num_symbols = transition_matrix.shape[0]
    seq_length = len(sequence)
    backward_arr = np.zeros((num_symbols, seq_length))
    # Starting probability will be considered equal for all states.
    start_prob = log(1 / 8.0)
    prev_symbols = NUC_TO_INDICES[sequence[-1]]
    for symbol in prev_symbols:
        backward_arr[symbol, -1] = start_prob
    # Now calculate for remainder of the sequence.
    for i in range(seq_length - 2, -1, -1):
        symbols = NUC_TO_INDICES[sequence[i]]
        for symbol in symbols:
            prob_transition_symbol1 = (
                backward_arr[prev_symbols[0], i+1] +
                transition_matrix[prev_symbols[0], symbol]
            )
            prob_transition_symbol2 = (
                backward_arr[prev_symbols[1], i+1] +
                transition_matrix[prev_symbols[1], symbol]
            )
            backward_arr[symbol, i] = (
                prob_transition_symbol1 +
                log(1 + exp(prob_transition_symbol2 -
                            prob_transition_symbol1))
            )
        prev_symbols = symbols

    final_transition_prob1 = (
        backward_arr[prev_symbols[0], 0] + start_prob)
    final_transition_prob2 = (
        backward_arr[prev_symbols[1], 0] + start_prob)
    final_probability = (
        final_transition_prob1 +
        log(1 + exp(final_transition_prob2 - final_transition_prob1))
    )
    return backward_arr, final_probability


def decode_posterior_probabilities(final_probability,
                                   forward_probabilities,
                                   backward_probabilities):
    """Calculates the posterior probability for each nucleotide in the
    sequence of being part of a CpG island.

    :param final_probability: overall probability of a sequence given
        the model paths
    :param forward_probabilities: pre-calculated forward probabilities
    :param backward_probabilities: pre-calculated backward probabilities
    :returns: a 1-dimensional array with the posterior probability of
        each nucleotide being in a CpG island

    """
    # We only care about the values in the first 4 rows, which represent
    # the CpG symbols.
    posterior_probabilities = ((forward_probabilities[:4] +
                                backward_probabilities[:4]).sum(0) -
                               final_probability)
    return np.exp(posterior_probabilities)


def identify_cpg_sites(posterior_probabilities, threshold=0.95):
    """Identifies probable CpG sites from posterior probabilities.

    :param posterior_probabilities: @todo
    :param threshold: @todo
    :returns: @todo

    """
    significant_positions = np.where(posterior_probabilities >=
                                     threshold)[0].tolist()
    sites = []
    start = prev_position = significant_positions[0]
    for position in significant_positions[1:]:
        if prev_position != position - 1:
            sites.append((start, prev_position))
            start = position
        prev_position = position
    if prev_position != position - 1:
        sites.append((start, prev_position))
    return sites


def output_probabilities(outfileh, probabilities, sequence):
    outlines = ["{},{}\n".format(nuc, prob) for nuc, prob in
                zip(sequence, probabilities)]
    outfileh.writelines(outlines)


def plot_probabilities(outfile_name, posterior_probabilities,
                       known_sites=None):
    plt.plot(posterior_probabilities, 'k-')
    if known_sites is not None:
        for start, end in known_sites:
            plt.axvspan(start, end, facecolor='r', alpha=0.3)
    plt.title('CpG site probabilities')
    plt.xlabel('sequence position')
    plt.ylabel('CpG site probability')
    plt.savefig(outfile_name)


def make_cli_parser():
    """Creates the command-line interface.

    :returns: an :py:class:`argparse.ArgumentParser` instance

    """
    cli_parser = argparse.ArgumentParser(description=__doc__)
    cli_parser.add_argument(
            'fasta_file', type=argparse.FileType('r'),
            help="FASTA format file containing one sequence"
    )
    cli_parser.add_argument(
            'transitions_file', type=argparse.FileType('r'),
            help=("CSV format file containing the Markov state "
                  "transitons")
    )
    cli_parser.add_argument(
            '--known-sites', type=argparse.FileType('r'),
            help=("CSV format file where each line contains the "
                  "start and end positions of a known CpG site")
    )
    return cli_parser


def main(argv=None):
    cli_parser = make_cli_parser()
    args = cli_parser.parse_args(argv)
    transition_matrix = np.log(np.loadtxt(args.transitions_file,
                                          delimiter=','))
    LOGGER.info("Parsing fasta file.")
    sequence = _parse_fasta_file(args.fasta_file)
    LOGGER.info("Calculating forward probabilities.")
    forward_probabilities, fwd_overall_prob = alglib.calc_forward_probabilities(
            sequence, transition_matrix)
    LOGGER.info("Calculating backward probabilities.")
    backward_probabilities, bkwd_overall_prob = calc_backward_probabilities(
            sequence, transition_matrix)
    LOGGER.info("Decoding sequence.")
    posterior_probabilities = decode_posterior_probabilities(
            fwd_overall_prob,
            forward_probabilities,
            backward_probabilities,
    )
    LOGGER.info("Writing probabilities to cpg_probabilities.csv")
    with open('cpg_probabilities.csv', 'w') as probabilities_outfile:
        output_probabilities(probabilities_outfile,
                             posterior_probabilities, sequence)
    LOGGER.info("Identifying CpG sites.")
    cpg_sites = identify_cpg_sites(posterior_probabilities)
    LOGGER.info("Writing CpG sites to cpg_sites.csv")
    with open('cpg_sites.csv', 'w') as sites_outfile:
        gencpgdata.output_sites(sites_outfile, cpg_sites)

    if args.known_sites:
        LOGGER.info("Parsing known sites file.")
        known_sites = _parse_known_sites_file(args.known_sites)
    else:
        known_sites = None
    LOGGER.info("Plotting probabilities to cpg_probabilities.png.")
    plot_probabilities('cpg_probabilities.png', posterior_probabilities,
                       known_sites)

if __name__ == '__main__':
    main()


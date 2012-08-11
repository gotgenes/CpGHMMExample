#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# Copyright (c) 2012 Christopher D. Lasher
#
# This software is released under the MIT License. Please see
# LICENSE.txt for details.


"""Algorithmic speedups for CpGHMMExample"""


import numpy as np

from libc.math cimport exp, log
cimport numpy as np


D_DTYPE = np.float64
ctypedef np.float64_t D_DTYPE_t


cpdef D_DTYPE_t log_add(D_DTYPE_t log_value_x, D_DTYPE_t log_value_y):
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


cdef void nuc_to_indices(char nucleotide, int indices[]):
    """Select the right matrix row indices for the given nucleotide.

    :param nucleotide: the current nucleotide of the sequence
    :param indices: a pre-allocated array of integers of length 2

    """
    if nucleotide == 'A':
        indices[0] = 0
    elif nucleotide == 'C':
        indices[0] = 1
    elif nucleotide == 'G':
        indices[0] = 2
    else:
        indices[0] = 3
    indices[1] = indices[0] + 4


def calc_forward_probabilities(
        char *sequence,
        np.ndarray[D_DTYPE_t, ndim=2] transition_matrix
    ):
    """Calculate the forward probabilities for a given sequence and
    transition matrix.

    :param sequence: @todo
    :param transition_matrix: @todo
    :returns: @todo

    """
    cdef int num_symbols = transition_matrix.shape[0]
    cdef long long seq_length = len(sequence)
    cdef np.ndarray[D_DTYPE_t, ndim=2] forward_arr = np.zeros(
            (num_symbols, seq_length))
    # Starting probability will be considered equal for all states.
    cdef D_DTYPE_t start_prob = log(1 / 8.0)
    cdef int prev_symbols[2]
    cdef int symbols[2]
    nuc_to_indices(sequence[0], prev_symbols)
    cdef int j
    for j in range(2):
        forward_arr[prev_symbols[j], 0] = start_prob
    # Now calculate for remainder of the sequence.
    cdef unsigned long i
    cdef int symbol
    cdef D_DTYPE_t prob_transition_symbol1, prob_transition_symbol2
    for i in range(1, seq_length):
        nuc_to_indices(sequence[i], symbols)
        for j in range(2):
            symbol = symbols[j]
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
        prev_symbols[0] = symbols[0]
        prev_symbols[1] = symbols[1]

    cdef D_DTYPE_t final_transition_prob1 = (
        forward_arr[prev_symbols[0], -1] + start_prob)
    cdef D_DTYPE_t final_transition_prob2 = (
        forward_arr[prev_symbols[1], -1] + start_prob)
    cdef D_DTYPE_t final_probability = (
        final_transition_prob1 +
        log(1 + exp(final_transition_prob2 - final_transition_prob1))
    )
    return forward_arr, final_probability


def calc_backward_probabilities(
        char *sequence,
        np.ndarray[D_DTYPE_t, ndim=2] transition_matrix
    ):
    """Calculate the backward probabilities for a given sequence and
    transition matrix.

    :param sequence: @todo
    :param transition_matrix: @todo
    :returns: @todo

    """
    cdef int num_symbols = transition_matrix.shape[0]
    cdef unsigned long seq_length = len(sequence)
    cdef np.ndarray[D_DTYPE_t, ndim=2] backward_arr = np.zeros(
            (num_symbols, seq_length))
    # Starting probability will be considered equal for all states.
    cdef D_DTYPE_t start_prob = log(1 / 8.0)
    cdef int prev_symbols[2]
    cdef int symbols[2]
    nuc_to_indices(sequence[seq_length-1], prev_symbols)
    cdef int j
    for j in range(2):
        backward_arr[prev_symbols[j], seq_length-1] = start_prob
    # Now calculate for remainder of the sequence.
    cdef long long i
    i = seq_length - 1
    cdef int symbol
    cdef D_DTYPE_t prob_transition_symbol1, prob_transition_symbol2
    for i in range(seq_length-2, -1, -1):
        nuc_to_indices(sequence[i], symbols)
        for j in range(2):
            symbol = symbols[j]
            prob_transition_symbol1 = (
                backward_arr[prev_symbols[0], i+1] +
                transition_matrix[symbol, prev_symbols[0]]
            )
            prob_transition_symbol2 = (
                backward_arr[prev_symbols[1], i+1] +
                transition_matrix[symbol, prev_symbols[1]]
            )
            backward_arr[symbol, i] = (
                prob_transition_symbol1 +
                log(1 + exp(prob_transition_symbol2 -
                            prob_transition_symbol1))
            )
        prev_symbols[0] = symbols[0]
        prev_symbols[1] = symbols[1]

    cdef D_DTYPE_t final_transition_prob1 = (
        backward_arr[prev_symbols[0], 0] + start_prob)
    cdef D_DTYPE_t final_transition_prob2 = (
        backward_arr[prev_symbols[1], 0] + start_prob)
    cdef D_DTYPE_t final_probability = (
        final_transition_prob1 +
        log(1 + exp(final_transition_prob2 - final_transition_prob1))
    )
    return backward_arr, final_probability


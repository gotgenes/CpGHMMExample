#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# Copyright (c) 2012 Christopher D. Lasher
#
# This software is released under the MIT License. Please see
# LICENSE.txt for details.


"""Tests for findcpgsites.py"""


from math import exp, log
import unittest

import numpy as np

from CpGHMMExample import findcpgsites
from CpGHMMExample.gencpgdata import TRANSITION_PROBABILITIES


LOG_PROBS = np.log(TRANSITION_PROBABILITIES)


class TestLogAdd(unittest.TestCase):
    """Tests for log_add()"""

    def test_log_add(self):
        cases = (
            (1, 2),
            (0.5, 0.5),
            (50, 0.7)
        )
        error = 0.0000000001
        for x, y in cases:
            expected = log(x + y)
            result = findcpgsites.log_add(log(x), log(y))
            difference = abs(expected - result)
            self.assertTrue(difference < error)
            result = findcpgsites.log_add(log(y), log(x))
            difference = abs(expected - result)
            self.assertTrue(difference < error)


class TestCalcForwardProbabilities(unittest.TestCase):
    """Tests for calc_forward_probabilities()"""

    def test_calc_forward_probabilities(self):
        sequence = 'CGCA'
        expected_prob_arr = np.zeros((8,4))
        expected_prob_arr[1,0] = expected_prob_arr[5,0] = log(1.0 / 8)

        expected_prob_arr[2,1] = findcpgsites.log_add(
                expected_prob_arr[1,0] + LOG_PROBS[1,2],
                expected_prob_arr[5,0] + LOG_PROBS[5,2]
        )
        expected_prob_arr[6,1] = findcpgsites.log_add(
                expected_prob_arr[1,0] + LOG_PROBS[1,6],
                expected_prob_arr[5,0] + LOG_PROBS[5,6]
        )

        expected_prob_arr[1,2] = findcpgsites.log_add(
                expected_prob_arr[2,1] + LOG_PROBS[2,1],
                expected_prob_arr[6,1] + LOG_PROBS[6,1]
        )
        expected_prob_arr[5,2] = findcpgsites.log_add(
                expected_prob_arr[2,1] + LOG_PROBS[2,5],
                expected_prob_arr[6,1] + LOG_PROBS[6,5]
        )

        expected_prob_arr[0,3] = findcpgsites.log_add(
                expected_prob_arr[1,2] + LOG_PROBS[1,0],
                expected_prob_arr[5,2] + LOG_PROBS[5,0]
        )
        expected_prob_arr[4,3] = findcpgsites.log_add(
                expected_prob_arr[1,2] + LOG_PROBS[1,4],
                expected_prob_arr[5,2] + LOG_PROBS[5,4]
        )

        expected_final_prob = findcpgsites.log_add(
                expected_prob_arr[0,3] + log(1.0 / 8),
                expected_prob_arr[4,3] + log(1.0 / 8)
        )

        result_prob_arr, result_final_prob = (
            findcpgsites.calc_forward_probabilities(
                    sequence, LOG_PROBS)
        )

        self.assertTrue(
                np.array_equal(result_prob_arr, expected_prob_arr))
        self.assertEqual(result_final_prob, expected_final_prob)


class TestCalcBackwardProbabilities(unittest.TestCase):
    """Tests for calc_backward_probabilities()"""

    def test_calc_backward_probabilities(self):
        sequence = 'CGCA'
        expected_prob_arr = np.zeros((8,4))
        expected_prob_arr[0,3] = expected_prob_arr[4,3] = log(1.0 / 8)
        expected_prob_arr[1,2] = findcpgsites.log_add(
                expected_prob_arr[0,3] + LOG_PROBS[1,0],
                expected_prob_arr[4,3] + LOG_PROBS[1,4]
        )
        expected_prob_arr[5,2] = findcpgsites.log_add(
                expected_prob_arr[0,3] + LOG_PROBS[5,0],
                expected_prob_arr[4,3] + LOG_PROBS[5,4]
        )

        expected_prob_arr[2,1] = findcpgsites.log_add(
                expected_prob_arr[1,2] + LOG_PROBS[2,1],
                expected_prob_arr[5,2] + LOG_PROBS[2,5]
        )
        expected_prob_arr[6,1] = findcpgsites.log_add(
                expected_prob_arr[1,2] + LOG_PROBS[6,1],
                expected_prob_arr[5,2] + LOG_PROBS[6,5]
        )

        expected_prob_arr[1,0] = findcpgsites.log_add(
                expected_prob_arr[2,1] + LOG_PROBS[1,2],
                expected_prob_arr[6,1] + LOG_PROBS[1,6]
        )
        expected_prob_arr[5,0] = findcpgsites.log_add(
                expected_prob_arr[2,1] + LOG_PROBS[5,2],
                expected_prob_arr[6,1] + LOG_PROBS[5,6]
        )

        expected_final_prob = findcpgsites.log_add(
                expected_prob_arr[1,0] + log(1.0 / 8),
                expected_prob_arr[5,0] + log(1.0 / 8)
        )

        result_prob_arr, result_final_prob = (
            findcpgsites.calc_backward_probabilities(
                    sequence, LOG_PROBS)
        )

        self.assertTrue(
                np.array_equal(result_prob_arr, expected_prob_arr))
        self.assertEqual(result_final_prob, expected_final_prob)


if __name__ == '__main__':
    unittest.main()


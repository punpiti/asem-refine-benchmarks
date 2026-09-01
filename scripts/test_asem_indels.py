"""Regression tests for the insertion-aware structural consensus update."""

from __future__ import annotations

import unittest

from asem_core import CombinedAlignmentCounts, accumulate_alignment, build_theta
from common import LocalAlignment


def aln(ref_aligned: str, read_aligned: str) -> LocalAlignment:
    return LocalAlignment(100.0, 0, 0, ref_aligned, read_aligned)


def consensus(theta: str, alignments: list[LocalAlignment]) -> str:
    counts = CombinedAlignmentCounts.empty(len(theta))
    for alignment in alignments:
        accumulate_alignment(counts, alignment)
    return build_theta(theta, counts, w=0.1)


class IndelConsensusTests(unittest.TestCase):
    def test_majority_insertion_adds_a_reference_column(self) -> None:
        inserted = aln("AC-GT", "ACAGT")
        unchanged = aln("ACGT", "ACGT")
        self.assertEqual(consensus("ACGT", [inserted] * 3 + [unchanged]), "ACAGT")

    def test_minority_insertion_is_discarded(self) -> None:
        inserted = aln("AC-GT", "ACAGT")
        unchanged = aln("ACGT", "ACGT")
        self.assertEqual(consensus("ACGT", [inserted] + [unchanged] * 2), "ACGT")

    def test_exactly_half_gap_retains_insertion(self) -> None:
        inserted = aln("AC-GT", "ACAGT")
        unchanged = aln("ACGT", "ACGT")
        self.assertEqual(consensus("ACGT", [inserted, unchanged]), "ACAGT")

    def test_multibase_insertion_preserves_aligned_order(self) -> None:
        inserted = aln("AC--GT", "ACTTGT")
        unchanged = aln("ACGT", "ACGT")
        self.assertEqual(consensus("ACGT", [inserted] * 3 + [unchanged]), "ACTTGT")

    def test_gap_majority_deletes_reference_column(self) -> None:
        deleted = aln("ACGT", "AC-T")
        unchanged = aln("ACGT", "ACGT")
        self.assertEqual(consensus("ACGT", [deleted] * 2 + [unchanged]), "ACT")

    def test_gap_plurality_below_half_does_not_delete(self) -> None:
        alignments = [
            aln("ACGT", "AC-T"),
            aln("ACGT", "AC-T"),
            aln("ACGT", "ACAT"),
            aln("ACGT", "ACCT"),
            aln("ACGT", "ACTT"),
        ]
        self.assertEqual(consensus("ACGT", alignments), "ACAT")


if __name__ == "__main__":
    unittest.main()

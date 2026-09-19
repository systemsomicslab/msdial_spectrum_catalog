"""The annotation kind comes from MS-DIAL's evidence column, and from the old prefix when there is none.

MS-DIAL used to say what an annotation rested on by prefixing the compound name: "no MS2: " for a
precursor-only suggestion, "low score: " for a spectrum that failed the search. Since 2026-09-15 it
writes the compound name alone and states the evidence in an "Evidence source" column.

Reading a table exported after that change through the prefixes alone classifies EVERY row as
msms_matched -- a precursor-only suggestion recorded as a reference match with a product-ion
spectrum. That is the overstatement classify_metabolite_name exists to prevent, and it is what this
catalog did between the MS-DIAL change and this one.

The column also says more than the prefix could. "low score: " collapsed two findings MS-DIAL now
separates -- compared and fell short, against compared and explained nothing -- and it had no way to
say a name came from an in-silico tool at all.
"""

from __future__ import annotations

import unittest

from msdial_spectrum_catalog.ingest import classify_metabolite_name


class EvidenceSourceDecidesTheKind(unittest.TestCase):
    def test_a_compared_spectrum_that_carried_the_match(self) -> None:
        for source in ("ReferenceSpectrum", "RuleBased"):
            kind, candidate, named = classify_metabolite_name("Quercetin", source)
            self.assertEqual("msms_matched", kind, source)
            self.assertEqual("Quercetin", candidate)
            self.assertTrue(named)

    def test_the_two_failures_the_prefix_collapsed_are_now_distinct(self) -> None:
        """"low score: " meant both of these, and they are opposite findings.

        A comparison that partly agreed is evidence for the candidate at a lower level. One that
        explained nothing is evidence against it. Recording them as one kind was the overstatement.
        """
        fell_short, _, _ = classify_metabolite_name("Quercetin", "WeakSpectrumMatch")
        explained_nothing, _, _ = classify_metabolite_name("Quercetin", "UnmatchedSpectrum")

        self.assertEqual("low_score", fell_short)
        self.assertEqual("unmatched_spectrum", explained_nothing)
        self.assertNotEqual(fell_short, explained_nothing)

    def test_a_precursor_only_suggestion_is_not_an_msms_match(self) -> None:
        """THE REGRESSION THIS FIXES. Without the column, this row read as msms_matched."""
        kind, candidate, _ = classify_metabolite_name("PC 34:1", "PrecursorOnly")

        self.assertEqual("precursor_only", kind)
        self.assertEqual("PC 34:1", candidate, "and the name is the compound name, with no prefix")

    def test_an_in_silico_name_is_neither_a_match_nor_a_bare_mass(self) -> None:
        """MS-FINDER, SIRIUS, CFM-ID and ICEBERG results, which the prefix could not express."""
        kind, _, _ = classify_metabolite_name("Quercetin", "InSilico")

        self.assertEqual("in_silico", kind)


class TheOldPrefixesStillWork(unittest.TestCase):
    """MS-DIAL 4 exports and MS-DIAL 5 exports made before the change stay ingestable."""

    def test_the_prefix_decides_when_no_column_was_exported(self) -> None:
        self.assertEqual("precursor_only", classify_metabolite_name("no MS2: PC 34:1")[0])
        self.assertEqual("precursor_only", classify_metabolite_name("w/o MS2: PEPTIDE")[0])
        self.assertEqual("low_score", classify_metabolite_name("low score: Quercetin")[0])
        self.assertEqual("msms_matched", classify_metabolite_name("Quercetin")[0])

    def test_the_prefix_is_still_stripped_from_the_candidate_name(self) -> None:
        _, candidate, _ = classify_metabolite_name("no MS2: PC 34:1")
        self.assertEqual("PC 34:1", candidate)

    def test_the_column_wins_over_a_prefix_if_a_file_somehow_carries_both(self) -> None:
        """Belt and braces: the column is what MS-DIAL measured, the prefix is a legacy rendering."""
        kind, candidate, _ = classify_metabolite_name("low score: Quercetin", "PrecursorOnly")

        self.assertEqual("precursor_only", kind)
        self.assertEqual("Quercetin", candidate)

    def test_an_unrecognised_or_absent_source_falls_through_rather_than_guessing(self) -> None:
        """Manual and Unspecified are deliberately unmapped: neither says what was compared.

        Falling through to the prefix means a table from either side of the MS-DIAL change is read
        by whatever evidence it does carry, and a term MS-DIAL adds later degrades to today's
        behaviour instead of vanishing into a wrong kind.
        """
        self.assertEqual("low_score", classify_metabolite_name("low score: X", "Manual")[0])
        self.assertEqual("low_score", classify_metabolite_name("low score: X", "null")[0])
        self.assertEqual("low_score", classify_metabolite_name("low score: X", "")[0])
        self.assertEqual("msms_matched", classify_metabolite_name("X", "SomethingNew")[0])


if __name__ == "__main__":  # pragma: no cover
    unittest.main()

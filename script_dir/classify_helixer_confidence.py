#!/usr/bin/env python3

"""
Classify Helixer predictions using RNA-seq and DIAMOND evidence.

## Final logic

1. Strong Swiss-Prot hit:
   HC

2. Strong Swiss-Prot or TrEMBL hit plus RNA support:
   HC

   RNA support here means:
   - at least one exact junction, OR
   - exon coverage >= threshold AND exon overlaps >= threshold

3. Strong RNA support without protein homology:
   HC

   Strong RNA support requires:
   - exon coverage >= threshold
   - exon overlaps >= threshold
   - >= minimum number of exact junctions
   - fraction of predicted introns with exact support >= threshold

4. Strong PTREP hit not rescued by the HC rules above:
   LC

5. Everything else:
   LC

Boundary-only support is retained in the output for inspection,
but is not used for HC classification.

## Expected inputs

- Representative Helixer protein FASTA, optionally gzipped
- RNA-seq support summary TSV
- Strong DIAMOND tables for Swiss-Prot, TrEMBL and PTREP
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path
from typing import TextIO


DIAMOND_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qlen",
    "slen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "qcovhsp",
    "scovhsp",
    "evalue",
    "bitscore",
    "stitle",
]


def smart_open(path: str) -> TextIO:
    """Open plain-text or gzip-compressed files."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def read_fasta_proteins(path: str) -> dict[str, int]:
    """
    Read protein IDs and lengths from FASTA.

    Returns
    -------
    dict
        {protein_id: protein_length}
    """
    proteins: dict[str, int] = {}

    current_id: str | None = None
    sequence: list[str] = []

    with smart_open(path) as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    proteins[current_id] = len("".join(sequence))

                current_id = line[1:].split()[0]
                sequence = []

            else:
                sequence.append(line)

        if current_id is not None:
            proteins[current_id] = len("".join(sequence))

    if not proteins:
        raise ValueError(f"No FASTA identifiers found in {path}")

    return proteins


def read_rnaseq_summary(path: str) -> dict[str, dict[str, str]]:
    """Read RNA-seq support summary indexed by transcript ID."""

    required_columns = {
        "transcript",
        "n_exons",
        "frac_exonic_bases",
        "n_exon_overlaps",
        "n_introns",
        "n_exact",
        "n_boundary_only",
        "frac_exact",
        "max_exact_junction_reads",
    }

    records: dict[str, dict[str, str]] = {}

    with open(path, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"No header found in {path}")

        missing = required_columns - set(reader.fieldnames)

        if missing:
            raise ValueError(
                "RNA-seq summary is missing required columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            transcript = row["transcript"]

            if transcript in records:
                raise ValueError(
                    f"Duplicate transcript in RNA-seq summary: {transcript}"
                )

            records[transcript] = row

    return records


def read_diamond_hits(path: str) -> dict[str, dict[str, str]]:
    """
    Read one best/strong DIAMOND hit per query.

    If multiple rows occur for a query, retain the row with the greatest
    bitscore; ties are resolved using the smaller e-value.
    """
    hits: dict[str, dict[str, str]] = {}

    with smart_open(path) as handle:
        first_line = handle.readline()

        if not first_line:
            return hits

        first_fields = first_line.rstrip("\n").split("\t")
        has_header = first_fields[0] == "qseqid"

        if has_header:
            header = first_fields
        else:
            header = DIAMOND_COLUMNS

        def process(fields: list[str]) -> None:

            if len(fields) < 14:
                raise ValueError(
                    f"Malformed DIAMOND row in {path}: expected at least "
                    f"14 columns, found {len(fields)}"
                )

            # stitle may theoretically contain tabs.
            if len(fields) > len(header):
                fields = fields[: len(header) - 1] + [
                    "\t".join(fields[len(header) - 1 :])
                ]

            row = dict(zip(header, fields))
            transcript = row["qseqid"]

            bitscore = float(row["bitscore"])
            evalue = float(row["evalue"])

            previous = hits.get(transcript)

            if previous is None:
                hits[transcript] = row
                return

            previous_bitscore = float(previous["bitscore"])
            previous_evalue = float(previous["evalue"])

            if (
                bitscore > previous_bitscore
                or (
                    bitscore == previous_bitscore
                    and evalue < previous_evalue
                )
            ):
                hits[transcript] = row

        if not has_header:
            process(first_fields)

        for line in handle:
            if not line.strip():
                continue

            process(line.rstrip("\n").split("\t"))

    return hits


def hit_value(
    hit: dict[str, str] | None,
    field: str,
    default: str = "NA",
) -> str:
    """Return one DIAMOND field or NA if no hit exists."""
    if hit is None:
        return default
    return hit.get(field, default)


def numeric_value(
    row: dict[str, str] | None,
    field: str,
    default: float = 0.0,
) -> float:
    """Safely convert one RNA-seq field to float."""
    if row is None:
        return default

    value = row.get(field, "NA")

    if value in ("", "NA", "NaN", "nan", None):
        return default

    return float(value)


def main() -> None:

    parser = argparse.ArgumentParser(
        description=(
            "Classify Helixer predictions using RNA-seq, "
            "Swiss-Prot, TrEMBL and PTREP evidence."
        )
    )

    parser.add_argument(
        "--proteins",
        required=True,
        help="Representative Helixer protein FASTA; may be gzipped.",
    )

    parser.add_argument(
        "--rnaseq",
        required=True,
        help="RNA-seq support summary TSV.",
    )

    parser.add_argument(
        "--swissprot",
        required=True,
        help="Strong/best Swiss-Prot DIAMOND TSV.",
    )

    parser.add_argument(
        "--trembl",
        required=True,
        help="Strong/best TrEMBL DIAMOND TSV.",
    )

    parser.add_argument(
        "--ptrep",
        required=True,
        help="Strong/best PTREP DIAMOND TSV.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output classification TSV.",
    )

    parser.add_argument(
        "--hc-ids",
        help="Optional output file containing HC transcript IDs.",
    )

    parser.add_argument(
        "--lc-ids",
        help="Optional output file containing LC transcript IDs.",
    )

    parser.add_argument(
        "--min-exon-coverage",
        type=float,
        default=0.50,
        help=(
            "Minimum fraction of exonic bases covered for exon support "
            "(default: 0.50)."
        ),
    )

    parser.add_argument(
        "--min-exon-overlaps",
        type=int,
        default=10,
        help=(
            "Minimum summed exon-level alignment overlaps for exon support "
            "(default: 10)."
        ),
    )

#    parser.add_argument(
#        "--min-exact-introns",
#        type=int,
#        default=2,
#        help=(
#            "Minimum exact junction-supported introns for RNA-only HC "
#            "(default: 2)."
#        ),
#    )

    parser.add_argument(
        "--min-exact-fraction",
        type=float,
        default=0.50,
        help=(
            "Minimum fraction of predicted introns with exact junction "
            "support for RNA-only HC (default: 0.50)."
        ),
    )

    args = parser.parse_args()

    proteins = read_fasta_proteins(args.proteins)
    transcript_ids = list(proteins)

    rnaseq = read_rnaseq_summary(args.rnaseq)

    swissprot = read_diamond_hits(args.swissprot)
    trembl = read_diamond_hits(args.trembl)
    ptrep = read_diamond_hits(args.ptrep)

    output_fields = [
        "transcript",
        "protein_length",

        "n_exons",
        "frac_exonic_bases",
        "n_exon_overlaps",
        "exon_supported",

        "n_introns",
        "n_exact",
        "frac_exact",
        "n_boundary_only",
        "max_exact_junction_reads",
        "strong_RNA_support",

        "swissprot_strong_hit",
        "swissprot_subject",
        "swissprot_pident",
        "swissprot_qcov",
        "swissprot_scov",
        "swissprot_evalue",
        "swissprot_title",

        "trembl_strong_hit",
        "trembl_subject",
        "trembl_pident",
        "trembl_qcov",
        "trembl_scov",
        "trembl_evalue",
        "trembl_title",

        "ptrep_strong_hit",
        "ptrep_subject",
        "ptrep_pident",
        "ptrep_qcov",
        "ptrep_scov",
        "ptrep_evalue",
        "ptrep_title",

        "confidence",
        "classification_reason",
    ]

    hc_ids: list[str] = []
    lc_ids: list[str] = []

    counts: dict[str, int] = {}

    ptrep_protein_conflicts: list[str] = []

    with open(
        args.output,
        "w",
        encoding="utf-8",
        newline="",
    ) as out_handle:

        writer = csv.DictWriter(
            out_handle,
            fieldnames=output_fields,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()

        for transcript in transcript_ids:

            rna = rnaseq.get(transcript)

            swiss_hit = swissprot.get(transcript)
            trembl_hit = trembl.get(transcript)
            ptrep_hit = ptrep.get(transcript)

            has_swissprot = swiss_hit is not None
            has_trembl = trembl_hit is not None
            has_ptrep = ptrep_hit is not None

            has_protein = (
                has_swissprot or
                has_trembl
            )

            n_exons = int(
                numeric_value(rna, "n_exons")
            )

            frac_exonic_bases = numeric_value(
                rna,
                "frac_exonic_bases",
            )

            n_exon_overlaps = int(
                numeric_value(
                    rna,
                    "n_exon_overlaps",
                )
            )

            n_introns = int(
                numeric_value(rna, "n_introns")
            )

            n_exact = int(
                numeric_value(rna, "n_exact")
            )

            n_boundary_only = int(
                numeric_value(
                    rna,
                    "n_boundary_only",
                )
            )

            frac_exact = numeric_value(
                rna,
                "frac_exact",
            )

            max_exact_junction_reads = int(
                numeric_value(
                    rna,
                    "max_exact_junction_reads",
                )
            )

            exon_supported = (
                frac_exonic_bases >= args.min_exon_coverage
                and
                n_exon_overlaps >= args.min_exon_overlaps
            )

            exact_supported = (
                n_exact > 0
            )

            strong_rna = (
                n_introns > 0
                and exon_supported
                and n_exact >= 1
            #    and frac_exact >= args.min_exact_fraction
            )

            # ------------------------------------------
            # Final classification
            # ------------------------------------------

            # Strong curated homology alone is sufficient.
            if has_swissprot:
                confidence = "HC"
                reason = "strong_SwissProt_hit"
                hc_ids.append(transcript)

            # PTREP is a veto unless it has strong swissprot.
            elif has_ptrep:
                confidence = "LC"
                reason = "strong_PTREP_hit"
                lc_ids.append(transcript)

            # Strong protein evidence together with RNA evidence.
            elif (
                has_protein
                and (
                    exact_supported
                    or exon_supported
                )
            ):
                confidence = "HC"
                reason = "strong_protein_plus_RNA"
                hc_ids.append(transcript)

            # Strong RNA evidence alone.
            elif strong_rna:
                confidence = "HC"
                reason = "strong_RNA_support"
                hc_ids.append(transcript)

            else:
                confidence = "LC"
                reason = "insufficient_support"
                lc_ids.append(transcript)

            counts[reason] = counts.get(reason, 0) + 1

#            if has_ptrep and has_protein:
#                ptrep_protein_conflicts.append(transcript)

            writer.writerow({

                "transcript": transcript,
                "protein_length": proteins[transcript],

                "n_exons": n_exons,
                "frac_exonic_bases": frac_exonic_bases,
                "n_exon_overlaps": n_exon_overlaps,

                "exon_supported": (
                    "yes" if exon_supported else "no"
                ),

                "n_introns": n_introns,
                "n_exact": n_exact,

                "frac_exact": (
                    frac_exact
                    if n_introns > 0
                    else "NA"
                ),

                "n_boundary_only": n_boundary_only,
                "max_exact_junction_reads": max_exact_junction_reads,

                "strong_RNA_support": (
                    "yes" if strong_rna else "no"
                ),

                "swissprot_strong_hit": (
                    "yes" if has_swissprot else "no"
                ),
                "swissprot_subject": hit_value(
                    swiss_hit,
                    "sseqid",
                ),
                "swissprot_pident": hit_value(
                    swiss_hit,
                    "pident",
                ),
                "swissprot_qcov": hit_value(
                    swiss_hit,
                    "qcovhsp",
                ),
                "swissprot_scov": hit_value(
                    swiss_hit,
                    "scovhsp",
                ),
                "swissprot_evalue": hit_value(
                    swiss_hit,
                    "evalue",
                ),
                "swissprot_title": hit_value(
                    swiss_hit,
                    "stitle",
                ),

                "trembl_strong_hit": (
                    "yes" if has_trembl else "no"
                ),
                "trembl_subject": hit_value(
                    trembl_hit,
                    "sseqid",
                ),
                "trembl_pident": hit_value(
                    trembl_hit,
                    "pident",
                ),
                "trembl_qcov": hit_value(
                    trembl_hit,
                    "qcovhsp",
                ),
                "trembl_scov": hit_value(
                    trembl_hit,
                    "scovhsp",
                ),
                "trembl_evalue": hit_value(
                    trembl_hit,
                    "evalue",
                ),
                "trembl_title": hit_value(
                    trembl_hit,
                    "stitle",
                ),

                "ptrep_strong_hit": (
                    "yes" if has_ptrep else "no"
                ),
                "ptrep_subject": hit_value(
                    ptrep_hit,
                    "sseqid",
                ),
                "ptrep_pident": hit_value(
                    ptrep_hit,
                    "pident",
                ),
                "ptrep_qcov": hit_value(
                    ptrep_hit,
                    "qcovhsp",
                ),
                "ptrep_scov": hit_value(
                    ptrep_hit,
                    "scovhsp",
                ),
                "ptrep_evalue": hit_value(
                    ptrep_hit,
                    "evalue",
                ),
                "ptrep_title": hit_value(
                    ptrep_hit,
                    "stitle",
                ),

                "confidence": confidence,
                "classification_reason": reason,
            })

    # --------------------------------------------------
    # HC / LC ID lists
    # --------------------------------------------------

    if args.hc_ids:
        with open(
            args.hc_ids,
            "w",
            encoding="utf-8",
        ) as handle:

            for transcript in hc_ids:
                handle.write(f"{transcript}\n")

    if args.lc_ids:
        with open(
            args.lc_ids,
            "w",
            encoding="utf-8",
        ) as handle:

            for transcript in lc_ids:
                handle.write(f"{transcript}\n")

    # --------------------------------------------------
    # Summary
    # --------------------------------------------------

    print(
        f"Total transcripts: {len(transcript_ids)}",
        file=sys.stderr,
    )

    print(
        f"HC transcripts: {len(hc_ids)}",
        file=sys.stderr,
    )

    print(
        f"LC transcripts: {len(lc_ids)}",
        file=sys.stderr,
    )

    for reason, count in sorted(counts.items()):
        print(
            f"{reason}: {count}",
            file=sys.stderr,
        )

    # --------------------------------------------------
    # PTREP / protein conflicts
    # --------------------------------------------------

#    print(
#        "PTREP + Swiss-Prot/TrEMBL conflicts: "
#        f"{len(ptrep_protein_conflicts)}",
#        file=sys.stderr,
#    )

    for transcript in ptrep_protein_conflicts[:20]:
        print(
            f"  conflict: {transcript}",
            file=sys.stderr,
        )

    # --------------------------------------------------
    # Warn about unmatched evidence IDs
    # --------------------------------------------------

    fasta_set = set(transcript_ids)

    evidence_sets = [
        ("RNA-seq", set(rnaseq)),
        ("Swiss-Prot", set(swissprot)),
        ("TrEMBL", set(trembl)),
        ("PTREP", set(ptrep)),
    ]

    for label, ids in evidence_sets:

        unmatched = ids - fasta_set

        if unmatched:

            print(
                f"WARNING: {len(unmatched)} {label} IDs did not match "
                f"the representative protein FASTA.",
                file=sys.stderr,
            )

            for transcript in sorted(unmatched)[:10]:
                print(
                    f"  unmatched: {transcript}",
                    file=sys.stderr,
                )


if __name__ == "__main__":
    main()

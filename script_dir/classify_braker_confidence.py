#!/usr/bin/env python3

"""
Combine BRAKER structural-support classification with DIAMOND results.

Final logic
-----------
1. Strong PTREP hit:
       LC

2. At least one transcript passes structural HC criteria:
       HC

3. No transcript passes structural HC criteria, but strong
   Swiss-Prot or TrEMBL hit:
       HC_homology_rescued

4. No structural HC transcript and no strong plant homology:
       LC

Expected inputs
---------------
- Representative BRAKER protein FASTA, optionally gzipped
- Structural HC transcript table from the gene-support filtering script
- Strong DIAMOND tables for Swiss-Prot, TrEMBL and PTREP

DIAMOND files may contain a header or may be headerless.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys


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


def normalize_braker_id(identifier: str) -> str:
    """
    Convert BRAKER protein/transcript IDs to gene IDs.

    Examples
    --------
    anno1.g1.t1                    -> g1
    anno1.both_agree_g10284.t1    -> both_agree_g10284
    g5265.t1                      -> g5265
    g5265                         -> g5265
    """

    identifier = identifier.strip().split()[0]

    # Remove annotation prefix used in BRAKER protein FASTA files.
    identifier = re.sub(r"^anno\d+\.", "", identifier)

    # Remove transcript suffix such as .t1, .t2, etc.
    identifier = re.sub(r"\.t\d+$", "", identifier)

    return identifier

def read_gff_gene_ids(path: str) -> list[str]:
    """Read all gene IDs from a GFF3 file."""

    gene_ids: list[str] = []
    seen: set[str] = set()

    with open(path, encoding="utf-8") as handle:
        for line in handle:

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            if fields[2] != "gene":
                continue

            attributes = fields[8]

            match = re.search(r"(?:^|;)ID=([^;]+)", attributes)

            if match is None:
                raise ValueError(
                    f"Gene feature without ID attribute:\n{line.rstrip()}"
                )

            gene_id = normalize_braker_id(match.group(1))

            if gene_id not in seen:
                seen.add(gene_id)
                gene_ids.append(gene_id)

    if not gene_ids:
        raise ValueError(f"No gene features found in {path}")

    return gene_ids


def read_structural_hc(path: str) -> set[str]:
    """
    Read gene IDs from the structural HC table.

    The input contains transcript-level rows that passed the structural
    support criteria. A gene is considered structurally HC if at least
    one of its transcripts is present.

    Accepts comment lines and a non-comment header beginning with gene_id.
    """

    hc_genes: set[str] = set()

    with open(path, encoding="utf-8") as handle:
        for line in handle:

            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")

            if fields[0] == "gene_id":
                continue

            gene_id = normalize_braker_id(fields[0])
            hc_genes.add(gene_id)

    return hc_genes


def read_diamond_hits(path: str) -> dict[str, dict[str, str]]:
    """
    Read one best/strong DIAMOND hit per query.

    If multiple rows occur for a query, retain the row with the greatest
    bitscore. Ties are resolved using the smaller e-value.
    """

    hits: dict[str, dict[str, str]] = {}

    with open(path, encoding="utf-8") as handle:

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

            # stitle may theoretically contain tabs; join extras back together.
            if len(fields) > len(header):
                fields = fields[: len(header) - 1] + [
                    "\t".join(fields[len(header) - 1 :])
                ]

            row = dict(zip(header, fields))
            gene_id = normalize_braker_id(row["qseqid"])

            bitscore = float(row["bitscore"])
            evalue = float(row["evalue"])

            previous = hits.get(gene_id)

            if previous is None:
                hits[gene_id] = row
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
                hits[gene_id] = row

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
    """Return a value from a DIAMOND hit, or a default if no hit exists."""

    if hit is None:
        return default

    return hit.get(field, default)


def main() -> None:

    parser = argparse.ArgumentParser(
        description=(
            "Combine BRAKER structural HC support with "
            "Swiss-Prot, TrEMBL and PTREP DIAMOND hits."
        )
    )

    parser.add_argument(
        "--gff",
        required=True,
        help="Final BRAKER4 GFF3; gene features define the complete gene set.",
    )

    parser.add_argument(
        "--structural-hc",
        required=True,
        help="HC transcript table generated from gene_support.tsv.",
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
        help="Optional output file containing HC gene IDs.",
    )

    parser.add_argument(
        "--lc-ids",
        help="Optional output file containing LC gene IDs.",
    )

    args = parser.parse_args()

    gene_ids = read_gff_gene_ids(args.gff)
    structural_hc = read_structural_hc(args.structural_hc)

    swissprot = read_diamond_hits(args.swissprot)
    trembl = read_diamond_hits(args.trembl)
    ptrep = read_diamond_hits(args.ptrep)

    output_fields = [
        "gene_id",
        "structural_HC",
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

        for gene_id in gene_ids:

            swiss_hit = swissprot.get(gene_id)
            trembl_hit = trembl.get(gene_id)
            ptrep_hit = ptrep.get(gene_id)

            is_structural_hc = gene_id in structural_hc

            has_plant_homology = (
                swiss_hit is not None
                or trembl_hit is not None
            )

            has_ptrep = ptrep_hit is not None

            # Strong PTREP similarity acts as a TE-related veto.
            if has_ptrep:

                confidence = "LC"
                reason = "strong_PTREP_hit"
                lc_ids.append(gene_id)

            # At least one transcript passed the structural HC criteria.
            elif is_structural_hc:

                confidence = "HC"
                reason = "structural_support"
                hc_ids.append(gene_id)

            # Structurally weak genes can be rescued by strong plant homology.
            elif has_plant_homology:

                confidence = "HC"
                reason = "homology_rescued"
                hc_ids.append(gene_id)

            # No structural HC transcript and no strong plant homology.
            else:

                confidence = "LC"
                reason = "low_structural_support_no_strong_homology"
                lc_ids.append(gene_id)

            counts[reason] = counts.get(reason, 0) + 1

            writer.writerow(
                {
                    "gene_id": gene_id,

                    "structural_HC": (
                        "yes" if is_structural_hc else "no"
                    ),

                    "swissprot_strong_hit": (
                        "yes" if swiss_hit is not None else "no"
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
                        "yes" if trembl_hit is not None else "no"
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
                        "yes" if ptrep_hit is not None else "no"
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
                }
            )

    # Sanity checks.
    if len(hc_ids) + len(lc_ids) != len(gene_ids):
        raise RuntimeError(
            f"Classification error: {len(gene_ids)} total genes, "
            f"but {len(hc_ids)} HC + {len(lc_ids)} LC."
        )

    overlap = set(hc_ids) & set(lc_ids)

    if overlap:
        raise RuntimeError(
            f"Classification error: {len(overlap)} genes occur "
            f"in both HC and LC."
        )

    if args.hc_ids:
        with open(
            args.hc_ids,
            "w",
            encoding="utf-8",
        ) as handle:

            for gene_id in hc_ids:
                handle.write(f"{gene_id}\n")

    if args.lc_ids:
        with open(
            args.lc_ids,
            "w",
            encoding="utf-8",
        ) as handle:

            for gene_id in lc_ids:
                handle.write(f"{gene_id}\n")

    print(
        f"Total genes: {len(gene_ids)}",
        file=sys.stderr,
    )

    print(
        f"HC genes: {len(hc_ids)}",
        file=sys.stderr,
    )

    print(
        f"LC genes: {len(lc_ids)}",
        file=sys.stderr,
    )

    for reason, count in sorted(counts.items()):
        print(
            f"{reason}: {count}",
            file=sys.stderr,
        )

    # Warn about evidence IDs not present in the representative FASTA.
    gff_set = set(gene_ids)

    for label, ids in [
        ("structural HC", structural_hc),
        ("Swiss-Prot", set(swissprot)),
        ("TrEMBL", set(trembl)),
        ("PTREP", set(ptrep)),
    ]:

        unmatched = ids - gff_set

        if unmatched:
            print(
                f"WARNING: {len(unmatched)} {label} IDs did not match "
                f"the retained GFF.",
                file=sys.stderr,
            )

            for gene_id in sorted(unmatched)[:10]:
                print(
                    f"  unmatched: {gene_id}",
                    file=sys.stderr,
                )


if __name__ == "__main__":
    main()

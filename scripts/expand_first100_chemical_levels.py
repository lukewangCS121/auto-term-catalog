#!/usr/bin/env python3
"""Expand first-100 chemical rows into name, concentration, and optimum rows."""

from __future__ import annotations

import argparse
import csv
import re
from collections.abc import Iterable
from pathlib import Path
from typing import NamedTuple
from urllib.parse import quote


class Level(NamedTuple):
    display: str
    evidence: str


LevelMap = dict[tuple[str, str, str], dict[str, list[Level]]]


def level_annotations() -> LevelMap:
    """Return source-audited levels keyed by document, subject, and chemical."""
    levels: LevelMap = {}

    def add(
        doc: int,
        subjects: Iterable[str],
        chemical: str,
        concentrations: Iterable[tuple[str, str]] = (),
        optima: Iterable[tuple[str, str]] = (),
    ) -> None:
        for subject in subjects:
            key = (str(doc), subject, chemical)
            if key in levels:
                raise ValueError(f"Duplicate level annotation: {key}")
            levels[key] = {
                "concentration": [Level(*item) for item in concentrations],
                "optimum": [Level(*item) for item in optima],
            }

    add(3, ["Ax23T"], "NaCl", [("3-6% NaCl", "3-6%")], [("optimum 3-4.5% NaCl", "3-4.5%")])
    add(3, ["Ax23T"], "H2", [("10 kPa H2", "10 kPa"), ("160 kPa H2", "160 kPa")])
    add(6, ["JP12T"], "NaCl", [("0-0.5% NaCl", "0-0.5%")], [("best 0.25% (w/v) NaCl", "0.25%")])
    add(7, ["P3C3T"], "NaCl", [("0.5% (v/w) NaCl", "0.5%")], [("optimum 0.5% (v/w) NaCl", "0.5%")])
    add(7, ["MAC6T"], "NaCl", [("4.0% (v/w) NaCl", "4.0%")], [("optimum 4.0% (v/w) NaCl", "4.0%")])
    add(15, ["AW1-3T"], "NaCl", [("0.0-0.5% (w/v) NaCl", "0.0-0.5%")])
    add(15, ["AW1-7T"], "NaCl", [("0-0.5% NaCl", "0-0.5%")])
    add(19, ["S174ᵀ", "W118ᵀ"], "NaCl", [("2.0% (w/v) NaCl", "2.0%")], [("optimum 2.0% (w/v) NaCl", "2.0%")])
    add(21, ["BD586T", "BD613T", "BD626T"], "NaCl", [("up to 10% NaCl", "10%")])
    add(41, ["YIM 135249T", "YIM 135347"], "sodium chloride", [("0-9.0% sodium chloride (NaCl, w/v)", "0-9.0%")], [("optimum 0% sodium chloride (NaCl, w/v)", "0%")])
    add(45, ["WXL103", "WXL210T"], "NaCl", [("2.5-8% (w/v) NaCl", "2.5-8%")], [("optimum 3-4% NaCl", "3-4%")])
    add(46, ["LJ205T", "TR449"], "NaCl", [("0.5-1.0% (w/v) NaCl", "0.5-1.0%")], [("optimum 0.5-1.0% (w/v) NaCl", "0.5-1.0%")])
    add(46, ["ZJ450T", "ZJ454", "ZJ70T", "ZJ77"], "NaCl", [("1.0% (w/v) NaCl", "1.0%")], [("optimum 1.0% (w/v) NaCl", "1.0%")])
    add(47, ["G4-2T"], "NaCl", [("0-1% (w/v) NaCl", "0-1%")], [("optimum 0.5% NaCl", "0.5%")])
    add(57, ["MINF-07-Sa-05T"], "NaCl", [("up to 9% (w/v) NaCl", "9%")])
    add(58, ["FAM 1755T"], "NaCl", [("0-5% (w/v) NaCl", "0-5%")])
    add(62, ["CHS3-5T"], "NaCl", [("2.0-11.0% NaCl", "2.0-11.0%")], [("optimum 3.0% NaCl", "3.0%")])
    add(62, ["M-2T"], "NaCl", [("2.0-4.0% NaCl", "2.0-4.0%")], [("optimum 3.0% NaCl", "3.0%")])
    add(63, ["SyP6RT"], "NaCl", [("0.25% (w/v) NaCl", "0.25%")], [("optimum 0.25% (w/v) NaCl", "0.25%")])
    add(64, ["CF4.4T"], "NaCl", [("0-20% (w/v) NaCl", "0-20%")], [("optimum 5% (w/v) NaCl", "5%")])
    add(64, ["KK5.5T"], "NaCl", [("0-20% (w/v) NaCl", "0-20%")], [("optimum 0-2.5% (w/v) NaCl", "0-2.5%")])
    add(65, ["Mg75T"], "NaCl", [("0-3% (w/v) NaCl", "0% to 3%")], [("optimal 0-1% (w/v) NaCl", "0-1%")])
    add(66, ["TMP9T", "TMP25"], "NaCl", [("up to 4.5% NaCl", "4.5%")])
    add(66, ["G.S.17T"], "NaCl", [("up to 3.0% NaCl", "3.0%")])
    add(69, ["REN36T"], "sodium chloride (NaCl)", [("0-2% sodium chloride (NaCl)", "0% to 2%")], [("optimal 0% sodium chloride (NaCl)", "salt concentration of 0%")])
    add(71, ["MAC3T"], "sodium chloride", [("0-15.0% (w/v) sodium chloride (NaCl)", "0-15.0%")])
    add(71, ["MAC8T"], "sodium chloride", [("0-14.0% (w/v) sodium chloride (NaCl)", "0-14.0%")])
    add(75, ["HSL-7T"], "NaCl", [("0-1% (w/v) NaCl", "0-1%")], [("optimum 0.5% NaCl", "0.5%")])
    add(78, ["28AT"], "sodium chloride", [("5-50 g l⁻¹ sodium chloride (NaCl)", "5 to 50 g l⁻¹")], [("optimum 20 g l⁻¹ sodium chloride (NaCl)", "20 g l⁻¹")])
    add(79, ["SD5T"], "sodium chloride", [("1-10% (w/v) sodium chloride", "1-10%")], [("optimum 5% sodium chloride", "5%")])
    add(84, ["DFM-14T"], "NaCl", [("2% (w/v) NaCl", "2%")], [("optimum 2% (w/v) NaCl", "2%")])
    add(88, ["AD34T", "PAK95"], "NaCl", [("1.4-4.8 M NaCl", "1.4-4.8 M")], [("optimum 3.1 M NaCl", "3.1 M")])
    add(90, ["N2T"], "NaCl", [("13-20% (w/v) NaCl", "13-20 %")], [("optimum 15% NaCl", "15%")])
    return levels


def sentence_spans(text: str):
    start = 0
    for match in re.finditer(r"[.!?](?=\s|$)", text):
        end = match.end()
        yield start, text[start:end]
        start = end
        while start < len(text) and text[start].isspace():
            start += 1
    if start < len(text):
        yield start, text[start:]


def evidence(row: dict[str, str], text: str, level: Level, level_type: str):
    base = row["label"]
    chemical_terms = [base]
    if base.casefold() in {"nacl", "sodium chloride", "sodium chloride (nacl)"}:
        chemical_terms = ["NaCl", "sodium chloride (NaCl)", "sodium chloride"]

    candidates = []
    amount_pattern = re.escape(level.evidence)
    if level.evidence[:1].isdigit():
        amount_pattern = rf"(?<![\d.\-]){amount_pattern}(?!\d)"
    for sentence_start, sentence in sentence_spans(text):
        amounts = list(re.finditer(amount_pattern, sentence, re.IGNORECASE))
        chemicals = [
            match
            for term in chemical_terms
            for match in re.finditer(re.escape(term), sentence, re.IGNORECASE)
        ]
        if not amounts or not chemicals:
            continue
        subject_present = bool(
            re.search(re.escape(row["relationship_subject_label"]), sentence, re.IGNORECASE)
        )
        for amount in amounts:
            prior = sentence[max(0, amount.start() - 45) : amount.start()]
            optimum_near = bool(re.search(r"(?i)(?:optimum|optimal|best)", prior))
            chemical = min(chemicals, key=lambda item: abs(item.start() - amount.start()))
            candidates.append(
                (
                    0 if subject_present else 1,
                    0 if optimum_near == (level_type == "optimum") else 1,
                    abs(chemical.start() - amount.start()),
                    sentence_start,
                    amount.start(),
                    chemical.start(),
                    sentence,
                    amount,
                    chemical,
                )
            )
    if not candidates:
        raise ValueError(
            f"Cannot locate evidence for doc {row['doc']} "
            f"{row['relationship_subject_label']!r} {level.display!r}"
        )

    *_, sentence, amount, chemical = min(candidates, key=lambda item: item[:6])
    sentence_start = min(candidates, key=lambda item: item[:6])[3]
    absolute_spans = sorted(
        {
            (sentence_start + amount.start(), sentence_start + amount.end()),
            (sentence_start + chemical.start(), sentence_start + chemical.end()),
        }
    )
    highlighted = sentence
    for start, end in sorted(
        {(amount.start(), amount.end()), (chemical.start(), chemical.end())}, reverse=True
    ):
        highlighted = highlighted[:start] + "[[" + highlighted[start:end] + "]]" + highlighted[end:]
    return (
        "; ".join(f"{start}:{end}" for start, end in absolute_spans),
        highlighted,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--documents-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    with args.input.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        old_fields = reader.fieldnames or []
        original_rows = list(reader)
    insert_at = old_fields.index("chemical_relationship") + 1
    new_fields = old_fields[:insert_at] + [
        "chemical_level_type",
        "chemical_base_label",
        "chebi_label",
    ] + old_fields[insert_at:]

    annotations = level_annotations()
    salt_template = next(
        row for row in original_rows if row["grounded_id"] == "CHEBI:26710"
    )
    texts = {
        path.name: path.read_text(encoding="utf-8")
        for path in args.documents_dir.glob("*.txt")
    }

    matched = set()
    output_rows = []
    for row in original_rows:
        relationship = row["field"] == "chemical_utilization_object"
        row["chemical_level_type"] = "chemical_name" if relationship else ""
        row["chemical_base_label"] = row["label"] if relationship else ""
        row["chebi_label"] = row["label"] if row["grounded_id"].startswith("CHEBI:") else ""
        key = (row["doc"], row["relationship_subject_label"], row["label"])
        annotation = annotations.get(key)
        if annotation:
            matched.add(key)
            if row["label"].casefold() in {"nacl", "sodium chloride", "sodium chloride (nacl)"}:
                candidates = [item for item in row["grounded_ids"].split("|") if item]
                row["grounded_id"] = "CHEBI:26710"
                row["grounded_ids"] = "|".join(dict.fromkeys(["CHEBI:26710", *candidates]))
                for field in ("kg_category", "kg_edge_count", "kg_edge_evidence"):
                    row[field] = salt_template[field]
                row["match_type"] = row["match_type"] or "synonym"
            row["kg_name"] = row["label"]
            row["chebi_label"] = row["label"]
        output_rows.append(row)

        if not annotation:
            continue
        for level_type in ("concentration", "optimum"):
            for level in annotation[level_type]:
                variant = dict(row)
                variant["entity_id"] = f"AUTO:{quote(level.display)}"
                variant["label"] = level.display
                variant["chemical_level_type"] = level_type
                variant["chemical_base_label"] = row["label"]
                variant["chebi_label"] = level.display
                variant["kg_name"] = level.display
                variant["match_type"] = f"context_{level_type}"
                variant["original_spans"], variant["context"] = evidence(
                    row, texts[row["source_file"]], level, level_type
                )
                output_rows.append(variant)

    missing = set(annotations) - matched
    if missing:
        raise ValueError(f"{len(missing)} annotated rows were not found: {sorted(missing)}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=new_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Original rows: {len(original_rows)}")
    print(f"Expanded rows: {len(output_rows)}")
    print(f"Annotated chemical relationships: {len(matched)}")


if __name__ == "__main__":
    main()

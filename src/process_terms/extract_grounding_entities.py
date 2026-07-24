#!/usr/bin/env python3
"""Flatten OntoGPT chemical-utilization YAML into grounding candidates."""

from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import quote, unquote

import yaml


FIELD_KINDS = {
    "study_taxa": "taxon_candidate",
    "strains": "strain",
    "chemicals_mentioned": "chemical",
}
RELATION_FIELDS = {
    "chemical_utilizations": ("chemical_utilization_subject", "strain", "chemical_utilization_object", "chemical"),
    "strain_relationships": ("strain_relationship_subject", "strain", "strain_relationship_object", "taxon_candidate"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--documents-dir",
        type=Path,
        help="Optional input directory whose sorted filenames supply source file and PMID provenance.",
    )
    return parser.parse_args()


def values(value: Any) -> Iterable[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return (item for item in value if isinstance(item, str))
    if isinstance(value, str):
        return [value]
    return []


def raw_values(raw_output: str, field: str) -> list[str]:
    match = re.search(
        rf"(?m)^[ \t]*{re.escape(field)}:[ \t]*([^\r\n]*)$",
        raw_output,
    )
    if not match or not match.group(1):
        return []
    return [value.strip() for value in match.group(1).split(";") if value.strip()]


def main() -> None:
    args = parse_args()
    source_files = (
        sorted(args.documents_dir.glob("*.txt"), key=lambda path: os.fsencode(path.name))
        if args.documents_dir
        else []
    )
    rows: list[dict[str, str]] = []

    with args.input.open(encoding="utf-8") as stream:
        documents = list(yaml.safe_load_all(stream))

    for doc_number, document in enumerate(documents, start=1):
        if not isinstance(document, dict):
            continue
        extracted = document.get("extracted_object") or {}
        entities = {
            entity.get("id"): entity.get("label", "")
            for entity in document.get("named_entities") or []
            if isinstance(entity, dict) and entity.get("id")
        }
        ids_by_label = {label.casefold(): entity_id for entity_id, label in entities.items() if label}
        source_file = source_files[doc_number - 1].name if doc_number <= len(source_files) else ""
        pmid_match = re.search(r"-(\d+)-abstract\.txt$", source_file)
        pmid = pmid_match.group(1) if pmid_match else ""
        seen: set[tuple[str, str]] = set()

        def add(field: str, kind: str, entity_id: str) -> None:
            key = (field, entity_id)
            if key in seen:
                return
            seen.add(key)
            rows.append(
                {
                    "doc": str(doc_number),
                    "source_file": source_file,
                    "pmid": pmid,
                    "field": field,
                    "kind": kind,
                    "entity_id": entity_id,
                    "label": entities.get(entity_id, unquote(entity_id.removeprefix("AUTO:"))),
                }
            )

        for field, kind in FIELD_KINDS.items():
            for entity_id in values(extracted.get(field)):
                add(field, kind, entity_id)
            for label in raw_values(document.get("raw_completion_output") or "", field):
                entity_id = ids_by_label.get(label.casefold(), f"AUTO:{quote(label)}")
                add(field, kind, entity_id)

        for relation_field, (subject_field, subject_kind, object_field, object_kind) in RELATION_FIELDS.items():
            for relation in extracted.get(relation_field) or []:
                if not isinstance(relation, dict):
                    continue
                if isinstance(relation.get("subject"), str):
                    add(subject_field, subject_kind, relation["subject"])
                if isinstance(relation.get("object"), str):
                    add(object_field, object_kind, relation["object"])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["doc", "source_file", "pmid", "field", "kind", "entity_id", "label"]
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()

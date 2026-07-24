#!/usr/bin/env python3
"""Ground extracted labels against merged-KG nodes and edge context."""

from __future__ import annotations

import argparse
import csv
import re
import unicodedata
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


CHEMICAL_FIELDS = {"chemicals_mentioned", "chemical_utilization_object"}
TAXON_FIELDS = {"study_taxa", "strain_relationship_object"}
INDICATOR_FIELDS = {
    "chemicals_utilized": "chemical_utilization_object",
    "study_taxa": "study_taxa",
    "strains": "strains",
}
DEFAULT_NODES = Path("/Users/lukewang/Downloads/merged-kg/merged-kg_nodes.tsv")
DEFAULT_EDGES = Path("/Users/lukewang/Downloads/merged-kg/merged-kg_edges.tsv")


@dataclass(frozen=True)
class Candidate:
    identifier: str
    name: str
    category: str
    match_type: str


def normalize(value: str) -> str:
    value = unicodedata.normalize("NFKC", value).casefold().strip()
    return re.sub(r"\s+", " ", value)


def entity_type(row: dict[str, str]) -> str:
    if row.get("field") in CHEMICAL_FIELDS or row.get("kind") == "chemical":
        return "chemical"
    if row.get("field") in TAXON_FIELDS or row.get("kind") == "taxon_candidate":
        return "taxon"
    if row.get("kind") == "strain":
        return "strain"
    return "other"


def category_matches(entity_kind: str, category: str) -> bool:
    if entity_kind == "chemical":
        return any(token in category for token in ("Chemical", "SmallMolecule", "MolecularEntity"))
    if entity_kind == "taxon":
        return "OrganismTaxon" in category
    if entity_kind == "strain":
        return "OrganismTaxon" in category
    return True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--nodes", default=DEFAULT_NODES, type=Path)
    parser.add_argument("--edges", default=DEFAULT_EDGES, type=Path)
    parser.add_argument("--max-edge-evidence", default=5, type=int)
    return parser.parse_args()


def candidate_rank(candidate: Candidate, entity_kind: str, edge_count: int) -> tuple[int, int, int, str]:
    prefix_preference = (
        (entity_kind == "chemical" and candidate.identifier.startswith("CHEBI:"))
        or (entity_kind in {"taxon", "strain"} and candidate.identifier.startswith("NCBITaxon:"))
    )
    return (
        0 if candidate.match_type == "name" else 1,
        0 if prefix_preference else 1,
        -edge_count,
        candidate.identifier,
    )


def main() -> None:
    args = parse_args()
    with args.input.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))

    wanted = {normalize(row.get("label", "")) for row in rows if row.get("label")}
    matches: dict[str, list[Candidate]] = defaultdict(list)
    with args.nodes.open(newline="", encoding="utf-8") as stream:
        for node in csv.DictReader(stream, delimiter="\t"):
            name_key = normalize(node.get("name", ""))
            if name_key in wanted:
                matches[name_key].append(Candidate(node["id"], node.get("name", ""), node.get("category", ""), "name"))
            for synonym in (node.get("synonym") or "").split("|"):
                synonym_key = normalize(synonym)
                if synonym_key in wanted and synonym_key != name_key:
                    matches[synonym_key].append(Candidate(node["id"], node.get("name", ""), node.get("category", ""), "synonym"))

    candidate_ids = {candidate.identifier for candidates in matches.values() for candidate in candidates}
    edge_counts: dict[str, int] = defaultdict(int)
    edge_evidence: dict[str, list[str]] = defaultdict(list)
    with args.edges.open(newline="", encoding="utf-8") as stream:
        for edge in csv.DictReader(stream, delimiter="\t"):
            subject = edge.get("subject", "")
            obj = edge.get("object", "")
            predicate = edge.get("predicate", "")
            for identifier, direction, neighbor in ((subject, "out", obj), (obj, "in", subject)):
                if identifier not in candidate_ids:
                    continue
                edge_counts[identifier] += 1
                if len(edge_evidence[identifier]) < args.max_edge_evidence:
                    edge_evidence[identifier].append(f"{direction}:{predicate}:{neighbor}")

    output_fields = list(rows[0].keys()) if rows else ["doc", "field", "kind", "entity_id", "label"]
    output_fields += list(INDICATOR_FIELDS)
    output_fields += [
        "grounded_id",
        "grounded_ids",
        "kg_name",
        "kg_category",
        "match_type",
        "kg_edge_count",
        "kg_edge_evidence",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=output_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            kind = entity_type(row)
            candidates = [
                candidate
                for candidate in matches.get(normalize(row.get("label", "")), [])
                if category_matches(kind, candidate.category)
            ]
            candidates = list({candidate.identifier: candidate for candidate in candidates}.values())
            candidates.sort(key=lambda candidate: candidate_rank(candidate, kind, edge_counts[candidate.identifier]))
            best = candidates[0] if candidates else None
            writer.writerow(
                {
                    **row,
                    **{
                        indicator: "1" if row.get("field") == source_field else "0"
                        for indicator, source_field in INDICATOR_FIELDS.items()
                    },
                    "grounded_id": best.identifier if best else "",
                    "grounded_ids": "|".join(candidate.identifier for candidate in candidates),
                    "kg_name": best.name if best else "",
                    "kg_category": best.category if best else "",
                    "match_type": best.match_type if best else "",
                    "kg_edge_count": str(edge_counts[best.identifier]) if best else "0",
                    "kg_edge_evidence": "|".join(edge_evidence[best.identifier]) if best else "",
                }
            )


if __name__ == "__main__":
    main()

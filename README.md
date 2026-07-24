# auto-term-catalog
code for extracting AUTO terms from ontoGPT output

## Generate a merged-KG-grounded TSV

The repository includes a deterministic command that converts a
multi-document OntoGPT extraction YAML file into a grounded TSV. It uses
`merged-kg_nodes.tsv` for exact normalized name/synonym matches and
`merged-kg_edges.tsv` for deterministic candidate ranking and edge evidence.

Requirements:

- Python 3
- PyYAML
- The extraction YAML
- The directory of source abstracts, named like
  `00001-41779015-abstract.txt`
- `merged-kg_nodes.tsv` and `merged-kg_edges.tsv`

Run:

```bash
./scripts/generate_merged_kg_grounded_tsv.sh \
  --extraction outputs/chemical_utilization_ijsem_first10_cborg_gpt41mini_no_grounding_final_20260608_145533.yaml \
  --documents-dir tmp/ijsem_first10_abstracts_with_pmids \
  --nodes /path/to/merged-kg_nodes.tsv \
  --edges /path/to/merged-kg_edges.tsv \
  --output outputs/chemical_utilization_ijsem_first10_merged_kg_grounded.tsv
```

The command also writes an intermediate `.entities.tsv` beside the output.
Use `--entities-output PATH` to choose another location. The process fixes its
locale, timezone, and Python hash seed; sorts source filenames bytewise; uses
stable candidate ordering; and writes outputs atomically. Given byte-identical
input files and the same command options, it produces byte-identical TSVs.

The grounded TSV includes binary `chemicals_utilized`, `study_taxa`, and
`strains` columns in addition to the selected ID, all candidate IDs, match
type, KG category, edge count, and representative edge evidence.

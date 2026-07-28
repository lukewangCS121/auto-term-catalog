#!/usr/bin/env bash
#
# Deterministically convert a multi-document OntoGPT extraction YAML file into
# a merged-KG-grounded TSV. The same ordered inputs produce byte-identical TSVs.

set -euo pipefail

export LC_ALL=C
export PYTHONHASHSEED=0
export TZ=UTC

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

usage() {
  cat <<'EOF'
Usage:
  scripts/generate_merged_kg_grounded_tsv.sh \
    --extraction PATH \
    --documents-dir PATH \
    --nodes PATH \
    --edges PATH \
    --metpo PATH \
    --output PATH \
    [--entities-output PATH] \
    [--max-edge-evidence NUMBER]

Required inputs:
  --extraction     Multi-document OntoGPT extraction YAML.
  --documents-dir Directory of PMID-labelled source .txt files. Files are
                  associated with YAML documents in bytewise filename order.
  --nodes          merged-kg_nodes.tsv.
  --edges          merged-kg_edges.tsv.
  --metpo          METPO ontology in RDF/XML OWL format.
  --output         Destination grounded TSV.

Optional:
  --entities-output     Destination for the intermediate flattened entity TSV.
                        Defaults to OUTPUT with ".entities.tsv" appended.
  --max-edge-evidence  Maximum incident KG edges retained per grounded ID.
                        Defaults to 5.
EOF
}

EXTRACTION=""
DOCUMENTS_DIR=""
NODES=""
EDGES=""
METPO=""
OUTPUT=""
ENTITIES_OUTPUT=""
MAX_EDGE_EVIDENCE=5

while (($#)); do
  case "$1" in
    --extraction)
      EXTRACTION="${2:?missing value for --extraction}"
      shift 2
      ;;
    --documents-dir)
      DOCUMENTS_DIR="${2:?missing value for --documents-dir}"
      shift 2
      ;;
    --nodes)
      NODES="${2:?missing value for --nodes}"
      shift 2
      ;;
    --edges)
      EDGES="${2:?missing value for --edges}"
      shift 2
      ;;
    --metpo)
      METPO="${2:?missing value for --metpo}"
      shift 2
      ;;
    --output)
      OUTPUT="${2:?missing value for --output}"
      shift 2
      ;;
    --entities-output)
      ENTITIES_OUTPUT="${2:?missing value for --entities-output}"
      shift 2
      ;;
    --max-edge-evidence)
      MAX_EDGE_EVIDENCE="${2:?missing value for --max-edge-evidence}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

for required in EXTRACTION DOCUMENTS_DIR NODES EDGES METPO OUTPUT; do
  if [[ -z "${!required}" ]]; then
    echo "Missing required argument: ${required}" >&2
    usage >&2
    exit 2
  fi
done

[[ -f "${EXTRACTION}" ]] || { echo "Extraction file not found: ${EXTRACTION}" >&2; exit 1; }
[[ -d "${DOCUMENTS_DIR}" ]] || { echo "Documents directory not found: ${DOCUMENTS_DIR}" >&2; exit 1; }
[[ -f "${NODES}" ]] || { echo "Node file not found: ${NODES}" >&2; exit 1; }
[[ -f "${EDGES}" ]] || { echo "Edge file not found: ${EDGES}" >&2; exit 1; }
[[ -f "${METPO}" ]] || { echo "METPO file not found: ${METPO}" >&2; exit 1; }
[[ "${MAX_EDGE_EVIDENCE}" =~ ^[0-9]+$ ]] || {
  echo "--max-edge-evidence must be a non-negative integer" >&2
  exit 2
}

if [[ -z "${ENTITIES_OUTPUT}" ]]; then
  ENTITIES_OUTPUT="${OUTPUT%.tsv}.entities.tsv"
fi

mkdir -p -- "$(dirname -- "${OUTPUT}")" "$(dirname -- "${ENTITIES_OUTPUT}")"
TMP_ENTITIES="$(mktemp "${ENTITIES_OUTPUT}.tmp.XXXXXX")"
TMP_OUTPUT="$(mktemp "${OUTPUT}.tmp.XXXXXX")"
cleanup() {
  rm -f -- "${TMP_ENTITIES}" "${TMP_OUTPUT}"
}
trap cleanup EXIT

python3 "${REPO_ROOT}/src/process_terms/extract_grounding_entities.py" \
  --input "${EXTRACTION}" \
  --documents-dir "${DOCUMENTS_DIR}" \
  --output "${TMP_ENTITIES}"

python3 "${REPO_ROOT}/src/process_terms/ground_entities_merged_kg.py" \
  --input "${TMP_ENTITIES}" \
  --nodes "${NODES}" \
  --edges "${EDGES}" \
  --metpo "${METPO}" \
  --max-edge-evidence "${MAX_EDGE_EVIDENCE}" \
  --strict-relationships \
  --output "${TMP_OUTPUT}"

chmod 0644 "${TMP_ENTITIES}" "${TMP_OUTPUT}"
mv -f -- "${TMP_ENTITIES}" "${ENTITIES_OUTPUT}"
mv -f -- "${TMP_OUTPUT}" "${OUTPUT}"

trap - EXIT
echo "Wrote entities: ${ENTITIES_OUTPUT}"
echo "Wrote grounded TSV: ${OUTPUT}"

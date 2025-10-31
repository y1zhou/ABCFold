#!/usr/bin/env bash
set -euo pipefail
if [[ "${TRACE-0}" == "1" ]]; then
    set -o xtrace
fi

DOCSTRING="Run ABCFold2 Boltz and Chai folding pipeline.

usage: ./$(basename "${0}") [options]

required arguments:
  -i, --input-yaml        input ABCFold2 YAML configuration file

options:
  -h, --help              show this help message and exit

  I/O
  ==========================================
  -o, --output-dir        output directory to store results (default: <mktemp -d>). Note
                            that content in this directory may be overwritten.

  -r, --run-id            the name of this run (default: basename of input YAML without extension)

  Steps
  ==========================================
  --no-search-msa         if set, skip MSA search run folding directly (default: false).
                            Note that if this is set, you probably also need to provide
                            '--chai-template-m8' and '--chai-template-cif-dir' to specify
                            the template information for Chai folding.

  --postprocess           if set, run postprocessing after folding (default: false)

  --dry-run               if set, perform a dry run without executing commands (default: false)

  Directory and file configs
  ==========================================
  --abcfold-repo-dir      path to the ABCFold2 repository (default: inferred from script location)

  --boltz-model-dir       path to the Boltz model checkpoints (default: \$BOLTZ_CACHE > ${HOME}/.boltz/)

  --chai-model-dir        path to the Chai model checkpoints (default: \$CHAI_DOWNLOADS_DIR > <chai-package-root>/downloads/)

  --template-cache-dir    path to the template cache directory (default: ${HOME}/.cache/rcsb/)

  --ccd-lib-dir           path to the CCD library directory (default: \$BOLTZ_CACHE/mols/)

  --chai-template-m8      path to the Chai template hits file (default: inferred from MSA search step)

  --chai-template-cif-dir path to the Chai template CIF files directory (default: inferred from MSA search step)
"
script_path="$(realpath "${0}")"
script_dir="$(dirname "${script_path}")"
abcfold_repo_dir="$(realpath "${script_dir}/..")"

template_cache_dir="${HOME}/.cache/rcsb"

search_msa=1
postprocess=0
dry_run=0

while [[ $# -gt 0 ]]; do
    case $1 in
    -h | --help)
        echo "${DOCSTRING}"
        exit 0
        ;;
    -i | --input-yaml)
        abcfold_config="$(realpath "${2}")"
        shift # past argument
        shift # past value
        ;;
    -o | --output-dir)
        abcfold_out_dir="${2}"
        shift
        shift
        ;;
    -r | --run-id)
        abcfold_run_id="${2}"
        shift
        shift
        ;;
    --no-search-msa)
        search_msa=0
        shift
        ;;
    --postprocess)
        postprocess=1
        shift
        ;;
    --dry-run)
        dry_run=1
        shift
        ;;
    --abcfold-repo-dir)
        abcfold_repo_dir="$(realpath "${2}")"
        shift
        shift
        ;;
    --boltz-model-dir)
        BOLTZ_CACHE="$(realpath "${2}")"
        export BOLTZ_CACHE
        shift
        shift
        ;;
    --chai-model-dir)
        CHAI_DOWNLOADS_DIR="$(realpath "${2}")"
        export CHAI_DOWNLOADS_DIR
        shift
        shift
        ;;
    --template-cache-dir)
        template_cache_dir="${2}"
        shift
        shift
        ;;
    --ccd-lib-dir)
        ccd_lib_dir="${2}"
        shift
        shift
        ;;
    --chai-template-m8)
        chai_template_m8="$(realpath "${2}")"
        shift
        shift
        ;;
    --chai-template-cif-dir)
        chai_template_cif_dir="$(realpath "${2}")"
        shift
        shift
        ;;
    *)
        echo "Unknown option: $1"
        echo "${DOCSTRING}"
        exit 1
        ;;
    esac
done

# Make sure required arguments are provided
if [ -z "${abcfold_config:-}" ]; then
    echo 'Error: input YAML file is required via '-i'.'
    echo "${DOCSTRING}"
    exit 1
fi

if [ ! -f "${abcfold_config}" ]; then
    echo "Error: input YAML file '${abcfold_config}' does not exist."
    exit 1
fi

# Set default run ID if not provided
if [ -z "${abcfold_run_id:-}" ]; then
    # Handle both .yaml and .yml extensions
    abcfold_run_id="$(basename "${abcfold_config%.*}")"
fi

# Set default output directory if not provided
if [ -z "${abcfold_out_dir:-}" ]; then
    abcfold_out_dir="$(mktemp -d)"
fi
abcfold_out_dir="$(realpath "${abcfold_out_dir}")"
mkdir -p "${abcfold_out_dir}"

# Set default Chai template paths
if [ -z "${chai_template_m8:-}" ]; then
    chai_template_m8="${abcfold_out_dir}/msa/all_chain_templates.m8"
fi
if [ -z "${chai_template_cif_dir:-}" ]; then
    chai_template_cif_dir="${abcfold_out_dir}/msa/chai_templates_cif"
fi

# Set model checkpoint cache directories
if [ -z "${BOLTZ_CACHE:-}" ]; then
    BOLTZ_CACHE="$(realpath "${HOME}/.boltz")"
fi
BOLTZ_CACHE="$(realpath "${BOLTZ_CACHE}")"
export BOLTZ_CACHE
ccd_lib_dir="${BOLTZ_CACHE}/mols"

if [ -z "${CHAI_DOWNLOADS_DIR:-}" ]; then
    echo '[NOTE] CHAI_DOWNLOADS_DIR not set, defaulting to <chai-package-root>/downloads/.'
else
    CHAI_DOWNLOADS_DIR="$(realpath "${CHAI_DOWNLOADS_DIR}")"
    export CHAI_DOWNLOADS_DIR
fi

# Print out settings before running
echo '=========================================='
echo 'ABCFold2 run settings:'
echo "  ABCFold2 repo:     ${abcfold_repo_dir}"
echo "  Input YAML:        ${abcfold_config}"
echo "  Output directory:  ${abcfold_out_dir}"
echo "  Run ID:            ${abcfold_run_id}"
echo "  Search MSA:        ${search_msa}"
echo "  Postprocess:       ${postprocess}"
echo "  Template cache:    ${template_cache_dir}"
echo "  CCD library:       ${ccd_lib_dir}"
echo "  Chai template M8 : ${chai_template_m8}"
echo "  Chai template CIF: ${chai_template_cif_dir}"
echo '=========================================='

if [ "${dry_run}" -eq 1 ]; then
    rmdir "${abcfold_out_dir}"
    exit 0
fi

# uv run abcfold2 validate $abcfold_config
cd "${abcfold_repo_dir}" || exit 1 # Change to repo dir to run uv commands

abcfold_msa_config="${abcfold_config}"
if [ "${search_msa}" -eq 1 ]; then
    echo '[PROGRESS] Running MSA search...'

    uv run abcfold2 prepare msa "${abcfold_config}" -f \
        -o "${abcfold_out_dir}" \
        --template-cache-dir "${template_cache_dir}"

    abcfold_msa_config="${abcfold_out_dir}/${abcfold_run_id}.yaml"
fi

echo '[PROGRESS] Preparing Boltz inputs...'
uv run abcfold2 prepare boltz "${abcfold_msa_config}" -o "${abcfold_out_dir}"
echo '[PROGRESS] Preparing Chai inputs...'
uv run abcfold2 prepare chai "${abcfold_msa_config}" -o "${abcfold_out_dir}" \
    --ccd-lib-dir "${ccd_lib_dir}"

echo '[PROGRESS] Folding input with Boltz...'
uv run abcfold2 fold boltz "${abcfold_msa_config}" \
    -i "${abcfold_out_dir}/boltz_${abcfold_run_id}/${abcfold_run_id}.yaml" \
    -o "${abcfold_out_dir}/boltz_${abcfold_run_id}/"

echo '[PROGRESS] Folding input with Chai...'
uv run abcfold2 fold chai "${abcfold_msa_config}" \
    -i "${abcfold_out_dir}/chai_${abcfold_run_id}/${abcfold_run_id}.yaml" \
    -o "${abcfold_out_dir}/chai_${abcfold_run_id}/" \
    --template-hits-path "${chai_template_m8}" \
    --template-cif-dir "${chai_template_cif_dir}"

if [ "${postprocess}" -eq 1 ]; then
    echo '[PROGRESS] Postprocessing results...'
    uv run abcfold2 postprocess collect -o "${abcfold_out_dir}" \
        -b "${abcfold_out_dir}/boltz_${abcfold_run_id}/" \
        -c "${abcfold_out_dir}/chai_${abcfold_run_id}/"
else
    echo "[PROGRESS] Skipping postprocessing as per user request.
Run the following command to postprocess later:

uv run abcfold2 postprocess collect -o \"${abcfold_out_dir}\" \\
    -b \"${abcfold_out_dir}/boltz_${abcfold_run_id}/\" \\
    -c \"${abcfold_out_dir}/chai_${abcfold_run_id}/\"
"
fi

echo '[SUCCESS] ABCFold2 run completed.'

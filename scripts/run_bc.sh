#!/usr/bin/env bash
set -euo pipefail
if [[ "${TRACE-0}" == "1" ]]; then
    set -o xtrace
fi

blue='\e[1;34m'
green='\e[1;32m'
end_hl='\e[0m'

print_cmd() {
    local col1="$1"
    local col2="$2"
    local col1_width="${3:-20}"
    printf "${blue}%-${col1_width}s${end_hl} %s\n" "$col1" "$col2"
}

DOCSTRING="Run ABCFold2 Boltz and Chai folding pipeline.

${green}Usage:${end_hl} ./$(basename "${0}") [options]

${green}Required arguments:${end_hl}
    $(print_cmd '-i, --input-yaml' 'input ABCFold2 YAML configuration file')

${green}Options:${end_hl}
    $(print_cmd '-h, --help' 'show this help message and exit')

    ${green}I/O${end_hl}
    ==========================================
    $(print_cmd '-o, --output-dir' 'output directory to store results (default: <mktemp -d>). Note that content in this directory may be overwritten.')
    $(print_cmd '-r, --run-id' 'the name of this run (default: basename of input YAML without extension)')

    ${green}Steps${end_hl}
    ==========================================
    $(print_cmd '--no-search-msa' 'if set, skip MSA search run folding directly (default: false). Note that if this is set, you probably also need to provide '--chai-template-m8' and '--chai-template-cif-dir' to specify the template information for Chai folding.')
    $(print_cmd '--msa-chains' 'if provided, only search MSA for the specified comma-separated chain IDs (default: all chains in input YAML). Example: A,B')
    $(print_cmd '--prepare-only' 'if set, only prepare inputs without running folding (default: false)')
    $(print_cmd '--no-fold-boltz' 'if set, skip Boltz folding step (default: false)')
    $(print_cmd '--no-fold-chai' 'if set, skip Chai folding step (default: false)')
    $(print_cmd '--postprocess' 'if set, run postprocessing after folding (default: false)')
    $(print_cmd '--dry-run' 'if set, perform a dry run without executing commands (default: false)')

    ${green}Directory and file configs${end_hl}
    ==========================================
    $(print_cmd '--abcfold-repo-dir' 'path to the ABCFold2 repository (default: inferred from script location)')
    $(print_cmd '--boltz-model-dir' "path to the Boltz model checkpoints (default: \$BOLTZ_CACHE > ${HOME}/.boltz/)")
    $(print_cmd '--chai-model-dir' "path to the Chai model checkpoints (default: \$CHAI_DOWNLOADS_DIR > <chai-package-root>/downloads/)")
    $(print_cmd '--template-cache-dir' "path to the template cache directory (default: ${HOME}/.cache/rcsb/)")
    $(print_cmd '--ccd-lib-dir' "path to the CCD library directory (default: \$BOLTZ_CACHE/mols/)")
    $(print_cmd '--chai-template-m8' 'path to the Chai template hits file (default: inferred from MSA search step)')
    $(print_cmd '--chai-template-cif-dir' 'path to the Chai template CIF files directory (default: inferred from MSA search step)')
"
script_path="$(realpath "${0}")"
script_dir="$(dirname "${script_path}")"
abcfold_repo_dir="$(realpath "${script_dir}/..")"

template_cache_dir="${HOME}/.cache/rcsb"

search_msa=1
prepare_only=0
skip_boltz=0
skip_chai=0
postprocess=0
dry_run=0

while [[ $# -gt 0 ]]; do
    case $1 in
    -h | --help)
        echo -e "${DOCSTRING}"
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
    --msa-chains)
        msa_chains="${2}"
        shift
        shift
        ;;
    --prepare-only)
        prepare_only=1
        shift
        ;;
    --no-fold-boltz)
        skip_boltz=1
        shift
        ;;
    --no-fold-chai)
        skip_chai=1
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
        echo -e "${DOCSTRING}"
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
conf_name="$(basename "${abcfold_config%.*}")"
if [ -z "${abcfold_run_id:-}" ]; then
    # Handle both .yaml and .yml extensions
    abcfold_run_id="${conf_name}"
fi

# Set default output directory if not provided
if [ -z "${abcfold_out_dir:-}" ]; then
    abcfold_out_dir="$(mktemp -d)"
fi
abcfold_out_dir="$(realpath "${abcfold_out_dir}")/${abcfold_run_id}"

# Set default Chai template paths
if [ -z "${chai_template_m8:-}" ]; then
    chai_template_m8="${abcfold_out_dir}/msa/all_chain_templates.m8"
fi
if [ -z "${chai_template_cif_dir:-}" ]; then
    chai_template_cif_dir="${abcfold_out_dir}/msa/templates"
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
print_cmd '  ABCFold2 repo:' "${abcfold_repo_dir}"
print_cmd '  Input YAML:' "${abcfold_config}"
print_cmd '  Output directory:' "${abcfold_out_dir}"
print_cmd '  Run ID:' "${abcfold_run_id}"
print_cmd '  Search MSA:' "${search_msa}"
print_cmd '  Prepare only:' "${prepare_only}"
print_cmd '  Skip Boltz fold:' "${skip_boltz}"
print_cmd '  Skip Chai fold:' "${skip_chai}"
print_cmd '  Postprocess:' "${postprocess}"
print_cmd '  Template cache:' "${template_cache_dir}"
print_cmd '  CCD library:' "${ccd_lib_dir}"
print_cmd '  Chai template M8:' "${chai_template_m8}"
print_cmd '  Chai template CIF:' "${chai_template_cif_dir}"
echo '=========================================='

if [ "${dry_run}" -eq 1 ]; then
    exit 0
fi

# uv run abcfold2 validate $abcfold_config
cd "${abcfold_repo_dir}" || exit 1 # Change to repo dir to run uv commands

mkdir -p "${abcfold_out_dir}"
abcfold_msa_config="${abcfold_out_dir}/${abcfold_run_id}.yaml"

if [ "${search_msa}" -eq 1 ]; then
    echo '[PROGRESS] Running MSA search...'

    abcfold_msa_args=('-f' '-o' "${abcfold_out_dir}" '--template-cache-dir' "${template_cache_dir}")
    if [ -n "${msa_chains:-}" ]; then
        abcfold_msa_args+=('--chains' "${msa_chains}")
    fi

    uv run abcfold2 prepare msa "${abcfold_config}" "${abcfold_msa_args[@]}"

    # in case $conf_name and $abcfold_run_id differ
    if [ "${conf_name}" != "${abcfold_run_id}" ]; then
        echo "[INFO] Moving MSA YAML to ${abcfold_msa_config}"
        mv "${abcfold_out_dir}/${conf_name}.yaml" "${abcfold_msa_config}"
    fi

    if [ ! -f "${abcfold_out_dir}/boltz_models/${abcfold_run_id}.yaml" ]; then
        echo '[PROGRESS] Preparing Boltz inputs...'
        uv run abcfold2 prepare boltz "${abcfold_msa_config}" -o "${abcfold_out_dir}"
    fi
    if [ ! -f "${abcfold_out_dir}/chai_models/${abcfold_run_id}.yaml" ]; then
        echo '[PROGRESS] Preparing Chai inputs...'
        uv run abcfold2 prepare chai "${abcfold_msa_config}" -o "${abcfold_out_dir}" \
            --ccd-lib-dir "${ccd_lib_dir}"
    fi
else
    echo '[PROGRESS] Skipping MSA search...'
    if [ ! -f "${abcfold_msa_config}" ]; then
        echo '[WARNING] YAML with MSA not found in output directory, copying input YAML as is.'
        cp -an "${abcfold_config}" "${abcfold_msa_config}"
    fi
fi

if [ "${prepare_only}" -eq 1 ]; then
    echo '[PROGRESS] Preparation only mode enabled, skipping folding steps.'
    exit 0
fi

if [ "${skip_boltz}" -eq 0 ]; then
    echo '[PROGRESS] Folding input with Boltz...'
    uv run abcfold2 fold boltz "${abcfold_msa_config}" \
        -i "${abcfold_out_dir}/boltz_models/${abcfold_run_id}.yaml" \
        -o "${abcfold_out_dir}/boltz_models/"
fi

if [ "${skip_chai}" -eq 0 ]; then
    echo '[PROGRESS] Folding input with Chai...'
    uv run abcfold2 fold chai "${abcfold_msa_config}" \
        -i "${abcfold_out_dir}/chai_models/${abcfold_run_id}.yaml" \
        -o "${abcfold_out_dir}/chai_models/" \
        --template-hits-path "${chai_template_m8}" \
        --template-cif-dir "${chai_template_cif_dir}"
fi

if [ "${postprocess}" -eq 1 ]; then
    echo '[PROGRESS] Postprocessing results...'
    uv run abcfold2 postprocess collect -o "${abcfold_out_dir}" \
        -b "${abcfold_out_dir}/boltz_models/" \
        -c "${abcfold_out_dir}/chai_models/"
else
    echo "[PROGRESS] Skipping postprocessing as per user request.
Run the following command to postprocess later:

uv run abcfold2 postprocess collect -o \"${abcfold_out_dir}\" \\
    -b \"${abcfold_out_dir}/boltz_models/\" \\
    -c \"${abcfold_out_dir}/chai_models/\"
"
fi

echo '[SUCCESS] ABCFold2 run completed.'

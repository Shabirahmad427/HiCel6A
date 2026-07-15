#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
src_root="$(cd "${repo_root}/.." && pwd)"
dest_root="${repo_root}/inputs"
manifest="${repo_root}/metadata/file_manifest.txt"

mkdir -p "${dest_root}" "$(dirname "${manifest}")"

if [ ! -f "${manifest}" ]; then
    printf 'Missing manifest: %s\n' "${manifest}" >&2
    exit 1
fi

while IFS= read -r archived_relpath; do
    [ -n "${archived_relpath}" ] || continue
    case "${archived_relpath}" in
        inputs/HiCel6A/*)
            src_relpath="1OCB/${archived_relpath#inputs/HiCel6A/}"
            dest_relpath="${archived_relpath#inputs/}"
            ;;
        inputs/HiCel6AD405N/*)
            src_relpath="1OC5/${archived_relpath#inputs/HiCel6AD405N/}"
            dest_relpath="${archived_relpath#inputs/}"
            ;;
        inputs/HiCel6AD416A/*)
            src_relpath="1OCJ/${archived_relpath#inputs/HiCel6AD416A/}"
            dest_relpath="${archived_relpath#inputs/}"
            ;;
        *)
            printf 'Unsupported manifest entry: %s\n' "${archived_relpath}" >&2
            exit 1
            ;;
    esac

    if [ ! -f "${src_root}/${src_relpath}" ]; then
        printf 'Missing source file for manifest entry: %s -> %s\n' "${archived_relpath}" "${src_relpath}" >&2
        exit 1
    fi

    mkdir -p "${dest_root}/$(dirname "${dest_relpath}")"
    cp -p "${src_root}/${src_relpath}" "${dest_root}/${dest_relpath}"
done < "${manifest}"

printf 'Collected %s files into %s\n' "$(wc -l < "${manifest}")" "${dest_root}"

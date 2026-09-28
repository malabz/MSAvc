#!/bin/bash
set -euo pipefail

cd "$SRC_DIR/src"

if [[ "${target_platform}" == osx-* ]]; then
    MAKE_FILE=Makefile.apple
else
    MAKE_FILE=Makefile
fi

make -f "$MAKE_FILE" -j"${CPU_COUNT:-1}" \
    CPP="${CXX}" \
    PREFIX="${PREFIX}"

mkdir -p "${PREFIX}/bin"
cp msavc_fasta vcf_merge "${PREFIX}/bin/"
cp msavc msavc_genome "${PREFIX}/bin/"
chmod +x "${PREFIX}/bin/msavc" "${PREFIX}/bin/msavc_genome"

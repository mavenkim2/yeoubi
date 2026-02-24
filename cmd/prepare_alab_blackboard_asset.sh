#!/usr/bin/env bash
set -euo pipefail

ENTRY_USD="${1:-/mnt/c/Users/maven/Downloads/ALab-2.2.0/ALab/entry.usda}"
BLACKBOARD_HINT="${2:-furniture_blackboard01_0001}"
OUT_DIR="${3:-assets/usd/blackboard}"
CAMERA_DISTANCE_SCALE="${4:-1.5}"

USD_BIN_DIR="src/third_party/install/usd/bin"
USDCAT="${USD_BIN_DIR}/usdcat"
ADD_CAMERA_BIN="build/yeoubi_usd_add_look_camera"

if [[ ! -x "${USDCAT}" ]]; then
  echo "Missing usdcat at ${USDCAT}" >&2
  exit 1
fi

if [[ ! -x "${ADD_CAMERA_BIN}" ]]; then
  echo "Missing camera tool at ${ADD_CAMERA_BIN}. Build first: cmake --build build --parallel" >&2
  exit 1
fi

mkdir -p "${OUT_DIR}"

RESOLVED_ENTRY_PRIM="$(${ADD_CAMERA_BIN} --resolve "${ENTRY_USD}" "${BLACKBOARD_HINT}")"
FLAT_USD="${OUT_DIR}/blackboard_flat.usda"
ASSET_USD="${OUT_DIR}/blackboard_asset.usda"

echo "Resolved blackboard prim: ${RESOLVED_ENTRY_PRIM}"
echo "Flatten+extract: ${ENTRY_USD}"
"${USDCAT}" --flatten --mask "${RESOLVED_ENTRY_PRIM}" "${ENTRY_USD}" -o "${FLAT_USD}"

RESOLVED_FLAT_PRIM="$(${ADD_CAMERA_BIN} --resolve "${FLAT_USD}" "${BLACKBOARD_HINT}")"
echo "Resolved extracted prim: ${RESOLVED_FLAT_PRIM}"

echo "Add camera"
"${ADD_CAMERA_BIN}" "${FLAT_USD}" "${RESOLVED_FLAT_PRIM}" "${ASSET_USD}" "${CAMERA_DISTANCE_SCALE}"

echo "Done"
echo "  flat:  ${FLAT_USD}"
echo "  asset: ${ASSET_USD}"

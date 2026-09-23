#!/bin/sh
# Generated from gen_all_multi.template.sh by setup.sh

# ---- Simulation name ----
SIMU_NAME="jupiter"

# -----------------------------------------------

MULTI_DIR="__MULTI_DIR__"
BIN_DIR="__BIN_DIR__"
TMP_DIR="__TMP_DIR__"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SIMU_DIR="${MULTI_DIR}/simu_${SIMU_NAME}"

# Create simulation directory.
mkdir -p "$SIMU_DIR"

# Generate paramsmulti.in from template and write into simulation directory.
PARAMS_TEMPLATE="${SCRIPT_DIR}/paramsmulti.template.in"
PARAMS_OUT="${SIMU_DIR}/paramsmulti.in"
if [ ! -f "$PARAMS_TEMPLATE" ]; then
    echo "ERROR: Template not found: $PARAMS_TEMPLATE" >&2
    exit 1
fi
sed -e "s|TMPDIR|${TMP_DIR}|g" \
    -e "s|SIMUNAME|${SIMU_NAME}|g" \
    "$PARAMS_TEMPLATE" > "$PARAMS_OUT"

# Move into simulation directory to run from there.
cd "$SIMU_DIR" || exit 1

# ---- Simulation parameters ----
${BIN_DIR}/gen_all_multi <<!
1   !   AU & Years
0   !   Invariable plane
0   !   No planet radii
1.   ! Mass of star
1   ! Number of planets
1.  ! Mass In Jupiter masses
5.2 0.03 0.000 0.00 0.00 0. ! a e i Omega omega M
plmulti.in
2.0
654876543
simu_${SIMU_NAME}
tpmulti.in
paramsmulti.in
100 ! # of particles
50    ! # of particles in first file
0. 0.01 ! Ex
3.  ! Max inc
2.4 2.6 !Amin & Amax
0
!
exit

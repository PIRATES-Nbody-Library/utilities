#!/bin/sh
# Generated from gen_all_multi.template.sh by setup.sh

# --- Set simulation name (set in paramsmulti and handed to gen_all_multi) ---
SIMU_NAME="jupiter"

# -----------------------------------------------

MULTI_DIR="__MULTI_DIR__"
BIN_DIR="__BIN_DIR__"
TMP_DIR="__TMP_DIR__"

# --- Simulation parameters (written to paramsmulti.in) ---
T0="0.0d0"
TSTOP="1.0e2"
DT="0.05"
DTOUT="1.0e1"
DTDUMP="1.0e1"
LFLAGS="F T F F T F"
RMIN="0.005"
RMAX="1000."
RMAXU="-1."
QMIN="-1."
LCLOSE="F"
FOPENSTAT="new"
# ============================================================

# Create simulation directory.
SIMU_DIR="${MULTI_DIR}/simu_${SIMU_NAME}"
mkdir -p "$SIMU_DIR"

# Generate paramsmulti.in directly in simulation directory.
cat > "${SIMU_DIR}/paramsmulti.in" <<PARAMS
${T0} ${TSTOP} ${DT}
${DTOUT} ${DTDUMP}
${LFLAGS}
${RMIN}  ${RMAX} ${RMAXU} ${QMIN} ${LCLOSE}
${TMP_DIR}
run_${SIMU_NAME}
bin_${SIMU_NAME}
${FOPENSTAT}
PARAMS

# Move into simulation directory to run from there.
cd "$SIMU_DIR" || exit 1

# --- Simulation parameters (directly piped to gen_all_multi) ---
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

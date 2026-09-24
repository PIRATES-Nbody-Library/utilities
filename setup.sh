#!/bin/sh
#
# setup.sh — Configuration and build script for PIRATES utilities
# Generate Fortran source, shell scripts, and makefile from templates
# with paths and compiler options, then optionally build executables.
#
# Usage:
#   ./setup.sh [options]
#
# Options:
#   --compiler=COMPILER        Fortran compiler (default: auto-detect)
#   --fflags=FLAGS             Compiler flags
#   --swift-dir=PATH           Path to swift_pirates/
#                              (default: ../Nbody_simulation_suite/swift_pirates)
#   --bin-dir=PATH             Executables directory (default: ./bin)
#   --tmp-dir=PATH             Temporary/dump dir (default: ./workspace/TMP)
#   --multi-dir=PATH           Multi-sim dir (default: ./workspace/multi)
#   --sbatch-account=NAME      SBATCH account (default: myaccount)
#   --sbatch-email=EMAIL       SBATCH email (default: user@example.com)
#   --sbatch-partition=NAME    SBATCH partition (default: standard)
#   --sbatch-time=TIME         SBATCH time limit (default: 240:00:00)
#   --sbatch-mem=MEM           SBATCH memory (default: 1gb)
#   --debug                    Enable debug flags
#   --build                    Build all utilities
#   --build=FILE               Build single file (name without .f)
#   --no-workspace             Skip workspace directory creation
#   --clean                    Remove generated .f files, executables, .o files
#   --help                     Show this help message
#

# ============================================================
# Default values
# ============================================================
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SWIFT_DIR=""
FORTRAN=""
FFLAGS=""
BIN_DIR=""
TMP_DIR=""
MULTI_DIR=""
SBATCH_ACCOUNT="myaccount"
SBATCH_EMAIL="user@example.com"
SBATCH_PARTITION="standard"
SBATCH_TIME="240:00:00"
SBATCH_MEM="1gb"
DO_BUILD=0
BUILD_TARGET=""
DO_CLEAN=0
DO_DEBUG=0
NO_WORKSPACE=0
OS=""

# ============================================================
# Parse command-line arguments
# ============================================================
print_help() {
    sed -n '/^# Usage:/,/^#$/p' "$0" | sed 's/^# \{0,1\}//'
    echo ""
    echo "Examples:"
    echo "  ./setup.sh                              # Auto-detect, generate files"
    echo "  ./setup.sh --compiler=gfortran          # Use gfortran"
    echo "  ./setup.sh --build                      # Generate files and build all"
    echo "  ./setup.sh --build=gen_all_multi         # Build single utility"
    echo "  ./setup.sh --debug --build              # Debug build"
    echo "  ./setup.sh --clean                      # Remove generated files"
    echo "  ./setup.sh --no-workspace               # Skip workspace creation"
}

for arg in "$@"; do
    case "$arg" in
        --compiler=*)
            FORTRAN="${arg#*=}"
            ;;
        --fflags=*)
            FFLAGS="${arg#*=}"
            ;;
        --swift-dir=*)
            SWIFT_DIR="${arg#*=}"
            ;;
        --bin-dir=*)
            BIN_DIR="${arg#*=}"
            ;;
        --tmp-dir=*)
            TMP_DIR="${arg#*=}"
            ;;
        --multi-dir=*)
            MULTI_DIR="${arg#*=}"
            ;;
        --sbatch-account=*)
            SBATCH_ACCOUNT="${arg#*=}"
            ;;
        --sbatch-email=*)
            SBATCH_EMAIL="${arg#*=}"
            ;;
        --sbatch-partition=*)
            SBATCH_PARTITION="${arg#*=}"
            ;;
        --sbatch-time=*)
            SBATCH_TIME="${arg#*=}"
            ;;
        --sbatch-mem=*)
            SBATCH_MEM="${arg#*=}"
            ;;
        --debug)
            DO_DEBUG=1
            ;;
        --build=*)
            DO_BUILD=1
            BUILD_TARGET="${arg#*=}"
            ;;
        --build)
            DO_BUILD=1
            ;;
        --no-workspace)
            NO_WORKSPACE=1
            ;;
        --clean)
            DO_CLEAN=1
            ;;
        --help|-h)
            print_help
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option: $arg"
            print_help
            exit 1
            ;;
    esac
done

# ============================================================
# Utility functions
# ============================================================
info()  { echo "[INFO]  $*"; }
warn()  { echo "[WARN]  $*"; }
error() { echo "[ERROR] $*" >&2; exit 1; }

# Resolve a path to absolute, creating it if needed.
resolve_path() {
    mkdir -p "$1" 2>/dev/null || error "Cannot create directory: $1"
    (cd "$1" && pwd)
}

# Resolve a path to absolute without creating it.
resolve_existing_path() {
    if [ ! -d "$1" ]; then
        error "Directory does not exist: $1"
    fi
    (cd "$1" && pwd)
}

# ============================================================
# Detect operating system
# ============================================================
detect_os() {
    case "$(uname -s)" in
        Linux*)  OS="linux"  ;;
        Darwin*) OS="macos"  ;;
        *)       OS="unknown" ;;
    esac
    info "Detected OS: $OS"
}

# ============================================================
# Detect Fortran compiler
# ============================================================
detect_compiler() {
    if [ -n "$FORTRAN" ]; then
        if ! command -v "$FORTRAN" > /dev/null 2>&1; then
            error "Specified compiler '$FORTRAN' not found in PATH"
        fi
        if ! "$FORTRAN" --version > /dev/null 2>&1; then
            error "Specified compiler '$FORTRAN' found but not functional"
        fi
        info "Using user-specified compiler: $FORTRAN"
        return
    fi

    for candidate in gfortran ifort ifx; do
        if command -v "$candidate" > /dev/null 2>&1; then
            if "$candidate" --version > /dev/null 2>&1; then
                FORTRAN="$candidate"
                info "Auto-detected compiler: $FORTRAN"
                return
            else
                warn "'$candidate' found in PATH but not functional — skipping"
            fi
        fi
    done

    error "No Fortran compiler found. Install gfortran or use --compiler=<path>"
}

# ============================================================
# Detect compiler version
# ============================================================
detect_compiler_version() {
    case "$FORTRAN" in
        gfortran*|ifort*|ifx*)
            COMPILER_VERSION=$("$FORTRAN" --version 2>&1 | head -1)
            ;;
        *)
            COMPILER_VERSION="unknown"
            ;;
    esac
    info "Compiler version: $COMPILER_VERSION"
}

# ============================================================
# Set compiler flags
# ============================================================
set_fflags() {
    if [ -n "$FFLAGS" ]; then
        info "Using user-specified FFLAGS: $FFLAGS"
        return
    fi

    if [ "$DO_DEBUG" -eq 1 ]; then
        case "$FORTRAN" in
            gfortran*)
                FFLAGS="-O0 -g -frecursive -ffixed-line-length-none -fbacktrace -ffpe-trap=zero,overflow,invalid -fcheck=all -Wall"
                ;;
            ifort*)
                FFLAGS="-O0 -g -recursive -extend-source -traceback -fpe0 -check all -warn all"
                ;;
            ifx*)
                FFLAGS="-O0 -g -recursive -extend-source -traceback -fpe0 -check all -warn all"
                ;;
            *)
                warn "Unknown compiler '$FORTRAN' — using generic debug flags"
                FFLAGS="-O0 -g"
                ;;
        esac
    else
        case "$FORTRAN" in
            gfortran*)
                FFLAGS="-O3 -frecursive -ffixed-line-length-none"
                ;;
            ifort*)
                FFLAGS="-O3 -recursive -extend-source"
                ;;
            ifx*)
                FFLAGS="-O3 -recursive -extend-source"
                ;;
            *)
                warn "Unknown compiler '$FORTRAN' — using generic flags"
                FFLAGS="-O3"
                ;;
        esac
    fi

    info "FFLAGS: $FFLAGS"
}

# ============================================================
# Resolve all paths
# ============================================================
resolve_paths() {
    # Swift directory — must exist and contain required files
    if [ -z "$SWIFT_DIR" ]; then
        SWIFT_DIR="${SCRIPT_DIR}/../Nbody_simulation_suite/swift_pirates"
    fi
    SWIFT_DIR="$(resolve_existing_path "$SWIFT_DIR")"

    # Validate swift directory contents
    if [ ! -f "$SWIFT_DIR/swift.inc" ]; then
        error "swift.inc not found in $SWIFT_DIR"
    fi
    if [ ! -f "$SWIFT_DIR/io/io.inc" ]; then
        error "io/io.inc not found in $SWIFT_DIR"
    fi
    if [ ! -f "$SWIFT_DIR/libswift.a" ]; then
        error "libswift.a not found in $SWIFT_DIR — build the simulation suite first"
    fi

    info "SWIFT_DIR: $SWIFT_DIR"

    # Bin directory
    if [ -z "$BIN_DIR" ]; then
        BIN_DIR="${SCRIPT_DIR}/bin"
    fi
    BIN_DIR="$(resolve_path "$BIN_DIR")"
    info "BIN_DIR: $BIN_DIR"

    # Workspace directories
    if [ "$NO_WORKSPACE" -eq 0 ]; then
        if [ -z "$TMP_DIR" ]; then
            TMP_DIR="${SCRIPT_DIR}/workspace/TMP"
        fi
        if [ -z "$MULTI_DIR" ]; then
            MULTI_DIR="${SCRIPT_DIR}/workspace/multi"
        fi

        for dir in "$TMP_DIR" "$MULTI_DIR"; do
            if [ -d "$dir" ]; then
                info "Directory already exists, skipping: $dir"
            fi
        done

        TMP_DIR="$(resolve_path "$TMP_DIR")"
        MULTI_DIR="$(resolve_path "$MULTI_DIR")"

        info "TMP_DIR: $TMP_DIR"
        info "MULTI_DIR: $MULTI_DIR"
    else
        # Still need values for template substitution
        if [ -z "$TMP_DIR" ]; then
            TMP_DIR="${SCRIPT_DIR}/workspace/TMP"
        fi
        if [ -z "$MULTI_DIR" ]; then
            MULTI_DIR="${SCRIPT_DIR}/workspace/multi"
        fi
        info "Skipping workspace directory creation (--no-workspace)"
    fi
}

# ============================================================
# Clean
# ============================================================
do_clean() {
    if [ "$DO_CLEAN" -eq 1 ]; then
        info "Cleaning..."

        # Remove generated Fortran source files (from templates)
        for f in gen_all_multi gen_multi_reprise; do
            target="$SCRIPT_DIR/fortran/${f}.f"
            [ -f "$target" ] && rm -f "$target" && info "Removed: $target"
        done

        # Remove object files
        find "$SCRIPT_DIR/fortran" -name "*.o" -delete 2>/dev/null

        # Remove generated makefile
        rm -f "$SCRIPT_DIR/fortran/@make_utilities"

        # Remove generated shell script
        rm -f "$SCRIPT_DIR/sh/gen_all_multi.sh"

        # Remove executables
        if [ -d "$SCRIPT_DIR/bin" ]; then
            rm -rf "$SCRIPT_DIR/bin"
            info "Removed bin/"
        fi

        # Remove macOS debug symbol bundles
        find "$SCRIPT_DIR" -name "*.dSYM" -type d \
            -exec rm -rf {} + 2>/dev/null

        info "Clean complete"
        exit 0
    fi
}

# ============================================================
# Generate Fortran source from templates
# ============================================================
generate_fortran() {
    info "Generating Fortran source files from templates..."

    for template in "$SCRIPT_DIR"/fortran/*.template.f; do
        [ -f "$template" ] || continue
        base="$(basename "$template" .template.f)"
        target="$SCRIPT_DIR/fortran/${base}.f"

        sed \
            -e "s|__SWIFT_DIR__|${SWIFT_DIR}|g" \
            -e "s|__BIN_DIR__|${BIN_DIR}|g" \
            -e "s|__TMP_DIR__|${TMP_DIR}|g" \
            -e "s|__MULTI_DIR__|${MULTI_DIR}|g" \
            -e "s|__SBATCH_ACCOUNT__|${SBATCH_ACCOUNT}|g" \
            -e "s|__SBATCH_EMAIL__|${SBATCH_EMAIL}|g" \
            -e "s|__SBATCH_PARTITION__|${SBATCH_PARTITION}|g" \
            -e "s|__SBATCH_TIME__|${SBATCH_TIME}|g" \
            -e "s|__SBATCH_MEM__|${SBATCH_MEM}|g" \
            "$template" > "$target"

        info "  ${base}.template.f -> ${base}.f"
    done
}

# ============================================================
# Generate shell scripts and input files from templates
# ============================================================
generate_shell() {
    info "Generating shell scripts and input files from templates..."

    # gen_all_multi.template.sh -> gen_all_multi.sh
    template="$SCRIPT_DIR/sh/gen_all_multi.template.sh"
    target="$SCRIPT_DIR/sh/gen_all_multi.sh"
    if [ -f "$template" ]; then
        sed \
        -e "s|__BIN_DIR__|${BIN_DIR}|g" \
        -e "s|__MULTI_DIR__|${MULTI_DIR}|g" \
        -e "s|__TMP_DIR__|${TMP_DIR}|g" \
        "$template" > "$target"
        chmod +x "$target"
        info "  gen_all_multi.template.sh -> gen_all_multi.sh"
    else
        warn "Template not found: $template"
    fi
}

# ============================================================
# Backup and generate @make_utilities
# ============================================================
generate_makefile() {
    template="$SCRIPT_DIR/fortran/@make_utilities.template"
    target="$SCRIPT_DIR/fortran/@make_utilities"

    if [ ! -f "$template" ]; then
        error "Template file not found: $template"
    fi

    # Backup existing @make_utilities
    if [ -f "$target" ]; then
        timestamp=$(date +%Y%m%d_%H%M%S)
        backup="${target}.backup.${timestamp}"
        cp "$target" "$backup"
        info "Backed up @make_utilities to $(basename "$backup")"
    fi

    date_str=$(date)

    sed \
        -e "s|__SWIFT_DIR__|${SWIFT_DIR}|g" \
        -e "s|__FORTRAN__|${FORTRAN}|g" \
        -e "s|__FFLAGS__|${FFLAGS}|g" \
        -e "s|__BIN_DIR__|${BIN_DIR}|g" \
        "$template" > "$target"

    # Replace template comment with generated-file header
    if [ "$OS" = "macos" ]; then
        sed -i '' \
            -e "s|^# @make_utilities template — used by setup.sh to generate @make_utilities|# Generated by setup.sh on ${date_str}|" \
            "$target"
    else
        sed -i \
            -e "s|^# @make_utilities template — used by setup.sh to generate @make_utilities|# Generated by setup.sh on ${date_str}|" \
            "$target"
    fi

    chmod +x "$target"
    info "Generated @make_utilities from template"
}

# ============================================================
# Build
# ============================================================
do_build() {
    if [ "$DO_BUILD" -eq 0 ]; then
        return
    fi

    if [ -n "$BUILD_TARGET" ]; then
        # Build a single file
        info "Building single target: $BUILD_TARGET"
        src="$SCRIPT_DIR/fortran/${BUILD_TARGET}.f"
        if [ ! -f "$src" ]; then
            error "Source file not found: $src"
        fi
        SWIFT_INC="$SWIFT_DIR"
        IO_INC="$SWIFT_DIR/io"
        echo "Compiling: $BUILD_TARGET"
        $FORTRAN $FFLAGS -I"$SWIFT_INC" -I"$IO_INC" \
            -o "$BIN_DIR/$BUILD_TARGET" "$src" \
            -L"$SWIFT_DIR" -lswift
        status=$?
        if [ $status -eq 0 ]; then
            info "Built: $BIN_DIR/$BUILD_TARGET"
        else
            error "Build failed for $BUILD_TARGET (exit code $status)"
        fi
    else
        # Build all
        info "Building all utilities..."
        cd "$SCRIPT_DIR/fortran" || error "Cannot cd to fortran/"
        sh "@make_utilities"
        status=$?
        if [ $status -eq 0 ]; then
            info "All utilities built successfully"
        else
            error "Build failed with exit code $status"
        fi
    fi
}

# ============================================================
# Print summary
# ============================================================
print_summary() {
    echo ""
    echo "========================================"
    echo "  PIRATES Utilities Configuration"
    echo "========================================"
    echo "  OS:              $OS"
    echo "  SWIFT_DIR:       $SWIFT_DIR"
    echo "  FORTRAN:         $FORTRAN"
    echo "  FFLAGS:          $FFLAGS"
    echo "  BIN_DIR:         $BIN_DIR"
    echo "  TMP_DIR:         $TMP_DIR"
    echo "  MULTI_DIR:       $MULTI_DIR"
    if [ "$DO_DEBUG" -eq 1 ]; then
        echo "  MODE:            debug"
    else
        echo "  MODE:            production"
    fi
    echo "========================================"
    echo ""
    if [ "$DO_BUILD" -eq 0 ]; then
        echo "To build all:      ./setup.sh --build"
        echo "To build one:      ./setup.sh --build=<name>"
        echo ""
        echo "Or manually:"
        echo "  cd fortran && sh @make_utilities"
    fi
}

# ============================================================
# Main
# ============================================================
main() {
    echo "========================================"
    echo "  setup.sh — Utilities Configuration"
    echo "========================================"
    echo ""

    detect_os

    # Handle clean early (before compiler detection)
    do_clean

    detect_compiler
    detect_compiler_version
    set_fflags
    resolve_paths

    generate_fortran
    generate_shell
    generate_makefile

    do_build

    print_summary
}

main "$@"

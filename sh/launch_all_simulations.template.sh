#!/bin/sh
# launch_all_simulations.template.sh -- Template for launching sub-simulations.
# Processed by gen_all_multi.sh.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

count=0
for script in go_rmvs3_*.sh; do
    [ -f "$script" ] || continue
    chmod +x "$script"
    echo "[INFO] Launching $script"
    ./"$script" > "${script%.sh}.log" 2>&1 &
    count=$((count + 1))
done

if [ "$count" -eq 0 ]; then
    echo "[ERROR] No go_rmvs3_*.sh scripts found." >&2
    exit 1
fi

echo "[INFO] Launched $count sub-simulation(s) in background."
echo "[INFO] Log files: go_rmvs3_*.log"

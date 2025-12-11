```sh
# Batch Solvation Tool *.pdb → solvated_*.pdb
# GROMACS can only be used on Linux systems
# Many proteins fail for unknown reasons
# Usage: Create an empty project in BASE_DIR on a Linux system with GROMACS installed, place this script in the project directory, and put all PDB files to be processed alongside it
# Remember to confirm that the gmx command is available first
# Fix: Remove dangerous grep cleanup, rely on GROMACS to automatically ignore HETATM

BASE_DIR="$HOME/gromacs_lab"
cd "$BASE_DIR"

LOG_DIR="$BASE_DIR/logs"
mkdir -p "$LOG_DIR"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$main_log"
}

main_log="$LOG_DIR/batch_$(date '+%Y%m%d_%H%M%S').log"
log "🚀 Starting batch solvation task... Working directory: $BASE_DIR"

# Capture errors but do not exit immediately (we control continue manually)
set +e

for pdb_file in [A-Za-z0-9][A-Za-z0-9][A-Za-z0-9][A-Za-z0-9].pdb; do
    [[ ! -f "$pdb_file" ]] && continue

    pdb_name="${pdb_file%.pdb}"
    output_pdb="solvated_${pdb_name}.pdb"

    if [[ -f "$output_pdb" ]]; then
        log "⏭️  $output_pdb already exists, skipping $pdb_file"
        continue
    fi

    log "🧪 Processing: $pdb_file → $output_pdb"

    work_dir="tmp_${pdb_name}_$$"
    mkdir -p "$work_dir"
    cp "$pdb_file" "$work_dir/"
    cd "$work_dir"

    task_log="$LOG_DIR/task_${pdb_name}.log"
    echo "=== Log started: $(date) ===" > "$task_log"

    # Modification: No manual grep! Use original PDB directly
    # GROMACS pdb2gmx will automatically ignore HETATM (ligands/water/ions etc.)

    log "   → pdb2gmx (force field 15, allow missing atoms)"
    if ! echo "15" | gmx pdb2gmx -f "${pdb_file}" -o protein.gro -water spce -missing >>"$task_log" 2>&1; then
        log "❌ Failed: $pdb_file failed in pdb2gmx step, check $task_log"
        cd "$BASE_DIR"
        rm -rf "$work_dir"
        continue
    fi

    log "   → editconf (build box)"
    if ! gmx editconf -f protein.gro -o box.gro -c -d 1.0 -bt cubic >>"$task_log" 2>&1; then
        log "❌ Failed: editconf error, check $task_log"
        cd "$BASE_DIR"
        rm -rf "$work_dir"
        continue
    fi

    log "   → solvate (add water)"
    if ! gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top >>"$task_log" 2>&1; then
        log "❌ Failed: solvate error, check $task_log"
        cd "$BASE_DIR"
        rm -rf "$work_dir"
        continue
    fi

    log "   → Convert to PDB"
    if ! gmx editconf -f solvated.gro -o "../${output_pdb}" >>"$task_log" 2>&1; then
        log "❌ Failed: PDB conversion error, check $task_log"
        cd "$BASE_DIR"
        rm -rf "$work_dir"
        continue
    fi

    cd "$BASE_DIR"
    rm -rf "$work_dir"
    log "✅ Completed: $output_pdb"

done

log "✨ Batch task completed!"

# Copy logs to Windows Downloads
windows_logs="/mnt/c/Users/sha1r/Downloads/logs" # Hardcoded the path because it's a hassle
if [[ -d "$windows_logs" ]]; then
    rm -rf "$windows_logs"
fi
cp -r "$LOG_DIR" "$windows_logs"
if [[ $? -eq 0 ]]; then
    log "✅ Logs copied to Windows Downloads: $windows_logs"
else
    log "❌ Failed to copy logs to Windows"
fi
```
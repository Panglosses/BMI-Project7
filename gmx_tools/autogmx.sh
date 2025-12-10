# 批量溶剂化小工具 *.pdb → solvated_*.pdb
# GROMACS只能在linux系统下使用
# 很多蛋白质不知为何跑不出来
# 使用方法：在安装好GROMACS的Linux系统的BASE_DIR处创建一个空项目，将改该脚本放进该项目目录中，并将所有待处理的pdb文件放在和它一起
# 记得先确认gmx指令是否已经可用
# 修复：移除危险的 grep 清理，依赖 GROMACS 自动忽略 HETATM

BASE_DIR="$HOME/gromacs_lab"
cd "$BASE_DIR"

LOG_DIR="$BASE_DIR/logs"
mkdir -p "$LOG_DIR"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$main_log"
}

main_log="$LOG_DIR/batch_$(date '+%Y%m%d_%H%M%S').log"
log "🚀 开始批量溶剂化任务... 工作目录: $BASE_DIR"

# 捕获错误但不立即退出（我们手动控制 continue）
set +e

for pdb_file in [A-Za-z0-9][A-Za-z0-9][A-Za-z0-9][A-Za-z0-9].pdb; do
    [[ ! -f "$pdb_file" ]] && continue

    pdb_name="${pdb_file%.pdb}"
    output_pdb="solvated_${pdb_name}.pdb"

    if [[ -f "$output_pdb" ]]; then
        log "⏭️  已存在 $output_pdb，跳过 $pdb_file"
        continue
    fi

    log "🧪 正在处理: $pdb_file → $output_pdb"

    work_dir="tmp_${pdb_name}_$$"
    mkdir -p "$work_dir"
    cp "$pdb_file" "$work_dir/"
    cd "$work_dir"

    task_log="$LOG_DIR/task_${pdb_name}.log"
    echo "=== 日志开始: $(date) ===" > "$task_log"

    # 修改：不再手动 grep！直接使用原始 PDB
    # GROMACS pdb2gmx 会自动忽略 HETATM（配体/水/离子等）

    log "   → pdb2gmx (力场 15, 允许缺失原子)"
    if ! echo "15" | gmx pdb2gmx -f "${pdb_file}" -o protein.gro -water spce -missing >>"$task_log" 2>&1; then
        log "❌ 失败: $pdb_file 在 pdb2gmx 步骤出错，请查看 $task_log"
        cd "$BASE_DIR"
        rm -rf "$work_dir"
        continue
    fi

    log "   → editconf (建盒子)"
    if ! gmx editconf -f protein.gro -o box.gro -c -d 1.0 -bt cubic >>"$task_log" 2>&1; then
        log "❌ 失败: editconf 出错，请查看 $task_log"
        cd "$BASE_DIR"
        rm -rf "$work_dir"
        continue
    fi

    log "   → solvate (加水)"
    if ! gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top >>"$task_log" 2>&1; then
        log "❌ 失败: solvate 出错，请查看 $task_log"
        cd "$BASE_DIR"
        rm -rf "$work_dir"
        continue
    fi

    log "   → 转换为 PDB"
    if ! gmx editconf -f solvated.gro -o "../${output_pdb}" >>"$task_log" 2>&1; then
        log "❌ 失败: 转换 PDB 出错，请查看 $task_log"
        cd "$BASE_DIR"
        rm -rf "$work_dir"
        continue
    fi

    cd "$BASE_DIR"
    rm -rf "$work_dir"
    log "✅ 完成: $output_pdb"

done

log "✨ 批量任务全部完成！"

# 复制 logs 到 Windows Downloads
windows_logs="/mnt/c/Users/sha1r/Downloads/logs" # 我嫌麻烦把路径硬编码了
if [[ -d "$windows_logs" ]]; then
    rm -rf "$windows_logs"
fi
cp -r "$LOG_DIR" "$windows_logs"
if [[ $? -eq 0 ]]; then
    log "✅ 日志已复制到 Windows Downloads: $windows_logs"
else
    log "❌ 复制日志到 Windows 失败"
fi

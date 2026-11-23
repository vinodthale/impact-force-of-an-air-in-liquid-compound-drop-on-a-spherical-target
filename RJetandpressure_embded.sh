#!/usr/bin/env bash
#SBATCH --job-name=force90
#SBATCH --partition=tyhcnormal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00

# In a case directory use: sbatch RJetandpressure_embded.sh
# On an allocated compute node: bash RJetandpressure_embded.sh
set -euo pipefail
export LC_ALL=C
case $# in
    0) ;;
    1) [[ $1 == --require-complete ]] || { echo 'Expected --require-complete or no arguments.' >&2; exit 2; } ;;
    *) echo 'Expected --require-complete or no arguments.' >&2; exit 2 ;;
esac

# A direct sbatch submission runs a copy of this script from the Slurm spool.
case_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
if [[ ! -f "$case_dir/Jetandpressure_embded.c" && -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    case_dir=$SLURM_SUBMIT_DIR
fi
cd -- "$case_dir"
case_dir=$(pwd -P)
db_dir=${case_dir##*/}
parent=${case_dir%/*}
vel_dir=${parent##*/}
if [[ ! "$vel_dir" =~ ^V[0-9]+\.[0-9]{2}$ || ! "$db_dir" =~ ^DB[0-9]+\.[0-9]{2}$ ]]; then
    echo "Run from a V<velocity>/DB<size> case directory: $case_dir" >&2
    exit 2
fi
for source in Jetandpressure_embded.c axi_embed_impact_force.h constants.h; do
    [[ -s "$source" ]] || { echo "Missing $case_dir/$source" >&2; exit 1; }
done
if [[ ! -d intermediate ]]; then
    echo "No saved intermediate/snapshot-* files in $case_dir" >&2
    exit 1
fi

work=$(mktemp -d -- "${TMPDIR:-/tmp}/axi_force_XXXXXXXX")
pending=''
cleanup() {
    if [[ -n "$pending" ]]; then rm -f -- "$pending"; fi
    rm -rf -- "$work"
}
trap cleanup EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM

find -L intermediate -maxdepth 1 -regextype posix-extended -type f \
    -regex 'intermediate/snapshot-[0-9]+\.[0-9]{4}' -print \
    | sort -t- -k2,2n > "$work/snapshots.list"
if [[ ! -s "$work/snapshots.list" ]]; then
    echo "No saved intermediate/snapshot-* files in $case_dir" >&2
    exit 1
fi

if command -v qcc >/dev/null && [[ -r "${BASILISK:-}/ast/interpreter/declarations.h" ]]; then
    basilisk_src=$(cd -- "$BASILISK" && pwd -P)
else
    basilisk_root=${BASILISK_ROOT:-$HOME/soft/basilisk-26.08.05-gcc12.2.0}
    basilisk_src="$basilisk_root/app/src"
    if [[ ! -r "$basilisk_root/scripts/env.sh" || ! -r "$basilisk_src/ast/interpreter/declarations.h" ]]; then
        echo "Set BASILISK to the configured Basilisk source directory or set BASILISK_ROOT to the cluster installation." >&2
        exit 1
    fi
    set +u
    module purge
    source "$basilisk_root/scripts/env.sh"
    module load compiler/openmpi/4.1.8-gcc12.2.0
    set -u
fi
ln -s -- "$basilisk_src" "$work/basilisk"
export BASILISK="$work/basilisk"
ln -s -- "$case_dir/intermediate" "$work/intermediate"
qcc -O2 -Wall -D_FORTIFY_SOURCE=0 Jetandpressure_embded.c -o "$work/extractor" -lm
(cd -- "$work" && ./extractor "${vel_dir#V}" "${db_dir#DB}" "$@")

# Publish the completed temporary file atomically on the case filesystem.
name="impact_force_Vel${vel_dir#V}_${db_dir}_90_deg.txt"
mkdir -p -- "$case_dir/Post"
pending=$(mktemp -- "$case_dir/Post/.impact_force_XXXXXXXX")
cp --preserve=mode -- "$work/$name" "$pending"
mv -fT -- "$pending" "$case_dir/Post/$name"
pending=''
printf '[OK] %s/Post/%s\n' "$case_dir" "$name"

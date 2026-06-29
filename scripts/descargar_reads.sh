#!/usr/bin/env bash
#
# descargar_reads.sh
# -----------------------------------------------------------------------------
# Descarga los reads (single-end) de los grupos 1-10 del practico de RNA-Seq
# desde el SRA/ENA usando sra-tools (prefetch + fasterq-dump).
#
# Todos los experimentos son SINGLE-END, por lo que fasterq-dump genera UN
# archivo .fastq por run (sin sufijos _1/_2).
#
# Uso:
#   ./descargar_reads.sh             # descarga TODOS los grupos (1 a 10)
#   ./descargar_reads.sh 7           # descarga solo el grupo 7
#   ./descargar_reads.sh 7 /ruta/destino
#
# Variables de entorno:
#   THREADS   numero de hebras para fasterq-dump   (default: 8)
#   DEST      carpeta destino de los .fastq         (default: $HOME/RNASEQ/reads)
#
# Requisitos: sra-tools (prefetch, fasterq-dump) en el PATH.
#   conda install -c bioconda sra-tools     # o el modulo del servidor
# -----------------------------------------------------------------------------
set -euo pipefail

THREADS="${THREADS:-8}"
DEST_DEFAULT="${DEST:-$HOME/RNASEQ/reads}"

# Grupo -> experimento (GXA/ArrayExpress), solo informativo
declare -A ACC=(
  [1]="E-GEOD-33294"   # carcinoma hepatocelular: tumor vs no-tumor adyacente
  [2]="E-GEOD-54505"   # MCF12A: miR-424 sobre-expresado vs vector vacio
  [3]="E-GEOD-62854"   # celulas estromales endometrio: PER2 knockdown vs control
  [4]="E-GEOD-58326"   # MCF7: siZNF217 vs control
  [5]="E-GEOD-54846"   # fibroblastos IMR90: shRNA macroH2A1 vs control
  [6]="E-GEOD-53280"   # HeLa + Salmonella Typhimurium
  [7]="E-MTAB-6013"    # iPSC vs cardiomiocito
  [8]="E-GEOD-42212"   # fibroblastos: senescencia inducida por Ras
  [9]="E-MTAB-8917"    # musculo esqueletico: diabetes T2 vs normal (4 vs 3)
  [10]="E-GEOD-52742"  # linfocitos B: asma alergica
)

# Grupo -> runs SINGLE-END (control ... test); 3 vs 3 salvo que se indique
declare -A RUNS=(
  [1]="SRR358994 SRR358996 SRR358998 SRR358995 SRR358997 SRR358999"
  [2]="SRR1146604 SRR1146605 SRR1146606 SRR1146607 SRR1146608 SRR1146609"
  [3]="SRR1635334 SRR1635335 SRR1635336 SRR1635337 SRR1635338 SRR1635339"
  [4]="SRR1363852 SRR1363853 SRR1363854 SRR1363855 SRR1363856 SRR1363857"
  [5]="SRR1166798 SRR1166799 SRR1166800 SRR1166801 SRR1166802 SRR1166803"
  [6]="SRR1049363 SRR1049364 SRR1049365 SRR1049366 SRR1049367 SRR1049368"
  [7]="ERR2365244 ERR2365245 ERR2365246 ERR2365247 ERR2365248 ERR2365249"
  [8]="SRR616151 SRR616152 SRR616153 SRR616154 SRR616155 SRR616156"
  [9]="ERR4010589 ERR4010590 ERR4010591 ERR4010592 ERR4010593 ERR4010594 ERR4010595"
  [10]="SRR1038580 SRR1038581 SRR1038582 SRR1038583 SRR1038584 SRR1038585"
)

# --- argumentos ---------------------------------------------------------------
GROUP="${1:-all}"
DEST="${2:-$DEST_DEFAULT}"

if [[ "$GROUP" == "all" ]]; then
  GROUPS=(1 2 3 4 5 6 7 8 9 10)
else
  if [[ -z "${RUNS[$GROUP]:-}" ]]; then
    echo "ERROR: grupo '$GROUP' no valido. Use uno de: 1 2 3 4 5 6 7 8 9 10 (o sin argumento para todos)." >&2
    exit 1
  fi
  GROUPS=("$GROUP")
fi

# --- chequeo de herramientas --------------------------------------------------
for tool in prefetch fasterq-dump; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "ERROR: no se encontro '$tool'. Instale sra-tools (conda install -c bioconda sra-tools)." >&2
    exit 1
  fi
done

mkdir -p "$DEST"
echo "Destino:  $DEST"
echo "Hebras:   $THREADS"
echo "Grupos:   ${GROUPS[*]}"
echo

# --- descarga -----------------------------------------------------------------
total=0; ok=0; skip=0; fail=0
for g in "${GROUPS[@]}"; do
  echo "=============================================================="
  echo " Grupo $g  ->  ${ACC[$g]}"
  echo "=============================================================="
  for run in ${RUNS[$g]}; do
    total=$((total+1))
    out="$DEST/${run}.fastq"
    if [[ -s "$out" || -s "${out}.gz" ]]; then
      echo "[$run] ya existe, se omite."
      skip=$((skip+1))
      continue
    fi
    echo "[$run] prefetch ..."
    if ! prefetch --max-size 100G -O "$DEST" "$run"; then
      echo "[$run] FALLO en prefetch." >&2
      fail=$((fail+1)); continue
    fi
    echo "[$run] fasterq-dump ..."
    if fasterq-dump --threads "$THREADS" --outdir "$DEST" "$DEST/$run/$run.sra" 2>/dev/null \
       || fasterq-dump --threads "$THREADS" --outdir "$DEST" "$run"; then
      rm -rf "$DEST/$run"   # limpia el .sra intermedio
      echo "[$run] OK -> ${run}.fastq"
      ok=$((ok+1))
    else
      echo "[$run] FALLO en fasterq-dump." >&2
      fail=$((fail+1))
    fi
  done
  echo
done

echo "=============================================================="
echo " Resumen: total=$total  ok=$ok  omitidos=$skip  fallidos=$fail"
echo "=============================================================="
[[ "$fail" -eq 0 ]]

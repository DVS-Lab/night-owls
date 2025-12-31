#!/usr/bin/env bash
set -euo pipefail

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

EVROOT="${maindir}/derivatives/fsl/EVFiles-FIR"

[[ -d "$EVROOT" ]] || { echo "ERROR: EV root not found: $EVROOT" >&2; exit 1; }

echo "Zeroing duration (2nd column) in EV txt files under:"
echo "  $EVROOT"

mapfile -d '' files < <(find "$EVROOT" -type f -name "*.txt" -print0)

if (( ${#files[@]} == 0 )); then
  echo "WARN: No .txt files found under: $EVROOT" >&2
  exit 0
fi

echo "Found ${#files[@]} .txt files."

changed=0
skipped=0

for f in "${files[@]}"; do
  bak="${f}.orig"
  [[ -e "$bak" ]] || cp -p "$f" "$bak"

  tmp="${f}.tmp"

  awk '
    /^[[:space:]]*$/ { print; next }
    /^[[:space:]]*#/ { print; next }
    {
      # Only modify lines that look like EV rows: >=3 fields
      if (NF >= 3) { $2 = 0 }
      print
    }
  ' "$bak" > "$tmp"

  mv "$tmp" "$f"
  changed=$((changed+1))
done

echo "Done. Updated: $changed files."
echo "Backups saved as: *.orig (next to each file)"

#!/usr/bin/env bash
# Link every skill in this repo into ~/.claude/skills/ so a git pull updates
# them in place. Symlinks, not copies -- the repo stays the single source of
# truth and the two can never drift.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$REPO/skills"
DEST="${CLAUDE_SKILLS_DIR:-$HOME/.claude/skills}"

uninstall=0
[[ "${1:-}" == "--uninstall" ]] && uninstall=1

[[ -d "$SRC" ]] || { echo "No skills/ directory in $REPO"; exit 1; }
mkdir -p "$DEST"

linked=0 skipped=0 removed=0

for skill in "$SRC"/*/; do
  [[ -d "$skill" ]] || continue
  name="$(basename "$skill")"
  target="$DEST/$name"

  if (( uninstall )); then
    # only remove links that point back into this repo -- never someone else's
    if [[ -L "$target" && "$(readlink "$target")" == "${skill%/}" ]]; then
      rm "$target"; echo "  removed  $name"; removed=$((removed+1))
    fi
    continue
  fi

  # An existing real directory is someone's own work -- never clobber it.
  if [[ -e "$target" && ! -L "$target" ]]; then
    echo "  SKIP     $name  (a real directory already exists at $target)"
    echo "           move or delete it first, then re-run"
    skipped=$((skipped+1)); continue
  fi

  [[ -L "$target" ]] && rm "$target"
  ln -s "${skill%/}" "$target"
  echo "  linked   $name"
  linked=$((linked+1))
done

echo
if (( uninstall )); then
  echo "Unlinked $removed skill(s) from $DEST"
else
  echo "Linked $linked skill(s) into $DEST${skipped:+, skipped $skipped}"
  echo "Restart Claude Code (or run /reload-plugins) to pick them up."
fi

# Remind about the data root, since every data-touching skill needs it.
if [[ -z "${LAB_DATA_ROOT:-}" ]] && (( ! uninstall )); then
  cat <<'EOF'

NOTE: LAB_DATA_ROOT is not set in this shell. Skills that read reference
data (genomes, GENCODE, GENIE, COSMIC) need it. Add to ~/.zshrc:

    export LAB_DATA_ROOT=/Users/Shared/lab-data/Database   # lab server
    export LAB_DATA_ROOT=~/Documents/GitHub/Database       # personal machine
EOF
fi

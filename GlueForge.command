#!/bin/sh
# GlueForge. Double-click this on macOS; a page opens in your browser and
# everything else happens there.
cd "$(dirname "$0")" || exit 1

if [ ! -f glueforge.py ]; then
  echo
  echo "  glueforge.py is not next to this file. Keep the folder together."
  echo
  read -r _ ; exit 1
fi

for py in python3 python; do
  if command -v "$py" >/dev/null 2>&1; then exec "$py" glueforge.py "$@"; fi
done

echo
echo "  Python 3 is not installed, or not on your PATH."
echo "  Get it from https://www.python.org/downloads/"
echo
read -r _

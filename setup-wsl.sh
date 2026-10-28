#!/usr/bin/env bash
# Set GlueForge up inside WSL (or any Debian/Ubuntu), including LightDock.
#
#   bash setup-wsl.sh            install everything, then start the server
#   bash setup-wsl.sh --check    say what is missing, change nothing
#   bash setup-wsl.sh --no-start install, do not start the server
#
# Why this file exists: LightDock publishes no Windows build for any Python
# version and its own page does not list Windows as supported. Its fastdfire
# scoring function is compiled C. On Linux that compiles in seconds once a
# toolchain is present, which is what this installs.
#
# The virtual environment goes in your Linux home directory, not next to the
# project. A venv on /mnt/c works but every import then crosses the Windows
# filesystem boundary, and on a package the size of scipy that is the
# difference between a second and most of a minute. The project itself stays
# exactly where it is; only the installed packages move.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV="${GLUEFORGE_VENV:-$HOME/.venvs/glueforge}"
APT_PACKAGES=(python3-venv python3-dev build-essential)

CHECK_ONLY=0
START=1
for arg in "$@"; do
  case "$arg" in
    --check)    CHECK_ONLY=1 ;;
    --no-start) START=0 ;;
    -h|--help)  sed -n '2,18p' "${BASH_SOURCE[0]}" | sed 's/^# \?//'; exit 0 ;;
    *) echo "unknown option: $arg" >&2; exit 2 ;;
  esac
done

say()  { printf '%s\n' "$*"; }
step() { printf '\n  [%s]\n' "$*"; }
ok()   { printf '    ok    %s\n' "$*"; }
bad()  { printf '    MISS  %s\n' "$*"; }

# --- refuse to run somewhere this cannot work ------------------------------

if [ "$(uname -s)" != "Linux" ]; then
  say "This script is for Linux or WSL. You are on $(uname -s)."
  say "On Windows, run GlueForge.bat instead -- but note that LightDock"
  say "cannot be installed there, so docking will be unavailable."
  exit 1
fi

if ! command -v apt-get >/dev/null 2>&1; then
  say "No apt-get here, so this script cannot install the toolchain."
  say "Install a C compiler, python3-venv and the Python headers the way"
  say "your distribution does it, then re-run with --no-start to continue."
  exit 1
fi

IN_WSL=0
grep -qiE 'microsoft|wsl' /proc/version 2>/dev/null && IN_WSL=1

say ""
say "  GlueForge setup"
say "  --------------------------------------------------------"
say "    project   $HERE"
say "    venv      $VENV"
say "    system    $( [ "$IN_WSL" = 1 ] && echo 'WSL' || echo 'Linux' )"

# --- what is already here ---------------------------------------------------

step "checking"
MISSING_APT=()
for p in "${APT_PACKAGES[@]}"; do
  if dpkg -s "$p" >/dev/null 2>&1; then ok "$p"; else bad "$p"; MISSING_APT+=("$p"); fi
done
command -v cc >/dev/null 2>&1 && ok "a C compiler" || bad "a C compiler"

HAVE_VENV=0
[ -x "$VENV/bin/python" ] && { ok "virtual environment"; HAVE_VENV=1; } \
                          || bad "virtual environment"

if [ "$HAVE_VENV" = 1 ]; then
  for m in numpy scipy lightdock; do
    if "$VENV/bin/python" -c "import $m" >/dev/null 2>&1; then ok "$m"; else bad "$m"; fi
  done
  if "$VENV/bin/python" -c "import rdkit" >/dev/null 2>&1; then
    ok "rdkit"
  else
    bad "rdkit (optional -- needed to grow chemistry)"
  fi
fi

if [ "$CHECK_ONLY" = 1 ]; then
  say ""
  say "  nothing was changed. Run without --check to install."
  exit 0
fi

# --- install ----------------------------------------------------------------

if [ "${#MISSING_APT[@]}" -gt 0 ]; then
  step "installing the toolchain"
  say "    sudo is needed once, for: ${MISSING_APT[*]}"
  sudo apt-get update -qq
  sudo apt-get install -y "${MISSING_APT[@]}"
fi

if [ "$HAVE_VENV" != 1 ]; then
  step "creating the virtual environment"
  python3 -m venv "$VENV"
fi

step "installing packages"
say "    lightdock compiles C here, so this takes a few minutes"
# Best-effort: an old pip still installs everything below, so a failure
# here is not worth stopping for, and under `set -e` it otherwise would.
"$VENV/bin/python" -m pip install --upgrade -q --retries 2 pip wheel \
  >/dev/null 2>&1 || true

"$VENV/bin/python" -m pip install --retries 2 numpy scipy || {
  say ""
  say "  Could not install numpy and scipy, which GlueForge needs for"
  say "  everything. These have builds for every platform, so this is"
  say "  almost certainly the network rather than the packages -- check"
  say "  that WSL has working DNS, and if you are behind a corporate"
  say "  proxy that pip has been told about it."
  exit 3
}
# LightDock separately, and without -q: it is the one that can fail, and when
# it does the output is the only thing that says why.
"$VENV/bin/python" -m pip install --retries 2 lightdock || {
  say ""
  say "  LightDock did not install. The output above says why; the two"
  say "  usual causes are a missing compiler (re-run this script, it"
  say "  installs build-essential) and a Python too new for one of its"
  say "  dependencies. Everything except docking works without it:"
  say ""
  say "      $VENV/bin/python $HERE/glueforge.py"
  exit 3
}
# Optional. A failure here costs the chemistry stage, not the docking stage,
# so it must not take the whole script down with it.
"$VENV/bin/python" -m pip install --retries 2 rdkit || say "    (rdkit failed; chemistry will be unavailable)"

# --- verify -----------------------------------------------------------------

step "verifying"
for t in lightdock3_setup.py lightdock3.py lgd_generate_conformations.py \
         lgd_cluster_bsas.py lgd_rank.py; do
  if [ -x "$VENV/bin/$t" ]; then ok "$t"; else bad "$t"; fi
done

say ""
PATH="$VENV/bin:$PATH" "$VENV/bin/python" "$HERE/sidecar/dock.py" --check || true

if [ "$START" = 0 ]; then
  say ""
  say "  done. To start GlueForge later:"
  say ""
  say "      PATH=\"$VENV/bin:\$PATH\" $VENV/bin/python $HERE/glueforge.py --no-browser"
  exit 0
fi

step "starting GlueForge"
if [ "$IN_WSL" = 1 ]; then
  say "    WSL forwards localhost, so the address below opens in your"
  say "    normal Windows browser. Ctrl-C here stops the server."
fi
say ""
# The venv's bin must be on PATH as well as being the interpreter: GlueForge
# shells out to the lightdock3*.py scripts by name.
exec env PATH="$VENV/bin:$PATH" "$VENV/bin/python" "$HERE/glueforge.py" --no-browser

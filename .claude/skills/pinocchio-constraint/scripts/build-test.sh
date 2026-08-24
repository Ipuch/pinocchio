#!/usr/bin/env bash
#
# Build (and optionally run) ONE pinocchio unit test out-of-tree, reusing the flags of an
# existing build directory. Much faster to iterate on than a full `cmake --build`, because
# it only ever compiles the single translation unit you are working on.
#
# Usage:
#   .claude/skills/pinocchio-constraint/scripts/build-test.sh unittest/my-constraint.cpp
#   ... --syntax-only     type-check only, no link, no run
#   ... --release         keep the build tree's -O3 -DNDEBUG instead of -UNDEBUG -O1 -g
#
# Env:
#   BUILD_DIR   build tree holding compile_commands.json (default: build)
#   OUT_DIR     where to put the binary                  (default: $TMPDIR or /tmp)
#   TEST_ARGS   extra arguments forwarded to the test binary
#
# By default asserts are ON (-UNDEBUG), which is what you want while developing a constraint:
# calc() guards are asserts. Re-check with --release before declaring victory.
#
set -euo pipefail

SRC=${1:?usage: build-test.sh <path/to/test.cpp> [--syntax-only] [--release]}
shift || true

SYNTAX_ONLY=0
RELEASE=0
for arg in "$@"; do
  case "$arg" in
    --syntax-only) SYNTAX_ONLY=1 ;;
    --release) RELEASE=1 ;;
    *) echo "unknown option: $arg" >&2; exit 2 ;;
  esac
done

ROOT=$(git -C "$(dirname "$SRC")" rev-parse --show-toplevel)
cd "$ROOT"

BUILD_DIR=${BUILD_DIR:-build}
DB="$BUILD_DIR/compile_commands.json"
[ -f "$DB" ] || {
  echo "no $DB — configure a build tree first (cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON)" >&2
  exit 1
}

OUT_DIR=${OUT_DIR:-${TMPDIR:-/tmp}}
NAME=$(basename "$SRC" .cpp)
BIN="$OUT_DIR/$NAME"
SER_DIR="$OUT_DIR/serialization-data"
mkdir -p "$OUT_DIR" "$SER_DIR"

# Reuse the compile line of any already-registered unittest as a flag template.
# NUL-separated argv, so tokens such as -DFOO="/a b/c" survive verbatim: do NOT switch this
# to `eval`, it eats the quotes that string-valued macros need.
mapfile -d '' -t CC < <(python3 - "$DB" <<'PY'
import json, shlex, sys

db = json.load(open(sys.argv[1]))
tpl = next((e for e in db if "/unittest/" in e["file"]), None) or db[0]
argv = shlex.split(tpl["command"])

out, skip_next = [], False
for a in argv[1:]:
    if skip_next:                       # the argument of -o
        skip_next = False
        continue
    if a == "-o":
        skip_next = True
        continue
    if a == "-c" or a.endswith(".cpp"):
        continue
    if a.startswith(("-DBOOST_TEST_MODULE=", "-DTEST_SERIALIZATION_FOLDER=")):
        continue                        # redefined by the caller below
    out.append(a)

sys.stdout.write("\0".join([argv[0]] + out) + "\0")
PY
)

CXX=${CC[0]}
FLAGS=("${CC[@]:1}"
       "-DBOOST_TEST_MODULE=${NAME//-/_}Test"
       "-DTEST_SERIALIZATION_FOLDER=\"$SER_DIR\"")
[ "$RELEASE" = 1 ] || FLAGS+=(-UNDEBUG -O1 -g)

if [ "$SYNTAX_ONLY" = 1 ]; then
  "$CXX" "${FLAGS[@]}" -fsyntax-only "$SRC"
  echo "syntax OK: $SRC"
  exit 0
fi

ENV_LIB="$(dirname "$(dirname "$CXX")")/lib"
"$CXX" "${FLAGS[@]}" "$SRC" -o "$BIN" \
  -L"$BUILD_DIR/src" -lpinocchio_default \
  -L"$ENV_LIB" -lboost_unit_test_framework -lboost_serialization -lboost_filesystem -lboost_system \
  -Wl,-rpath,"$ROOT/$BUILD_DIR/src" -Wl,-rpath,"$ENV_LIB"

echo "built: $BIN"
if [ -n "${TEST_ARGS:-}" ]; then
  # shellcheck disable=SC2086
  exec "$BIN" $TEST_ARGS
fi
exec "$BIN"

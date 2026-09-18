#!/bin/zsh
# A01: build the wheel, install into a FRESH environment, run from a path
# containing spaces with PYTHONPATH cleared and the research repository not
# importable. Fails if the old repo could ever satisfy the import.
set -e
TOOLKIT="$(cd "$(dirname "$0")/.." && pwd)"
FRESH="$(mktemp -d /tmp/tk_fresh_XXXX)"
WORK="$(mktemp -d '/tmp/tk work dir with spaces XXXX')"
cd "$TOOLKIT"
"${PYTHON:-python3}" -m venv "$FRESH/venv"
"$FRESH/venv/bin/pip" -q install --upgrade pip >/dev/null
"$FRESH/venv/bin/pip" -q build 2>/dev/null || "$FRESH/venv/bin/pip" -q install build
"$FRESH/venv/bin/python" -m build --wheel --outdir "$FRESH/dist" >/dev/null
WHEEL="$(ls "$FRESH"/dist/*.whl | head -1)"
echo "wheel: $(basename "$WHEEL")"
"$FRESH/venv/bin/pip" -q install "$WHEEL"

# copy ONLY toolkit examples into the spaced working directory
mkdir -p "$WORK/examples"
cp "$TOOLKIT"/examples/*.yaml "$TOOLKIT"/examples/*.csv "$WORK/examples/"

cd "$WORK"
env -u PYTHONPATH "$FRESH/venv/bin/python" - <<'EOF'
import tweezer_experiment, json, sys
loc = tweezer_experiment.__file__
assert "site-packages" in loc, f"imported from unexpected place: {loc}"
assert "monte-carlo" not in loc.lower() and "level2_joint_transfer" not in loc, loc
print("import location:", loc)
b = tweezer_experiment.ExperimentBundle.load("examples/short_transfer.yaml")
r = tweezer_experiment.simulate(b, shots=96)
print("simulate from spaced cwd: final S =", float(r.survival[-1]), "unbound =", r.initial_unbound)
m = tweezer_experiment.export_controls(b, "export_out")
print("export + readback:", json.dumps(m["readback"]["physical_csv_roundtrip"]))
sys.exit(0)
EOF
echo "A01 OK"

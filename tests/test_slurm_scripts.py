import glob
import os
import re
import subprocess

import pytest
from conftest import ROOT

SLURM_DIR = os.path.join(ROOT, "slurm")
SCRIPTS = sorted(glob.glob(os.path.join(SLURM_DIR, "*.sh")))

# The active Hellbender scripts (everything except the historical Longleaf one).
HELLBENDER = [s for s in SCRIPTS if "longleaf" not in os.path.basename(s)]


def test_scripts_present():
    names = {os.path.basename(s) for s in SCRIPTS}
    assert {
        "glide-dock.sh",
        "ligprep.sh",
        "glide-array.sh",
        "prepwizard.sh",
        "grid-gen.sh",
        "longleaf-glide-gpu.sh",
    } <= names


@pytest.mark.parametrize("script", SCRIPTS, ids=[os.path.basename(s) for s in SCRIPTS])
def test_bash_syntax(script):
    # Every script must at least pass `bash -n`.
    r = subprocess.run(["bash", "-n", script], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


@pytest.mark.parametrize("script", HELLBENDER, ids=[os.path.basename(s) for s in HELLBENDER])
def test_hellbender_real_values(script):
    text = open(script).read()
    assert "--account dkhf3-lab" in text
    assert "--partition dkhf3-lab" in text
    assert "module load schrodinger/v2024-2" in text


def test_array_script_uses_array_task_id():
    text = open(os.path.join(SLURM_DIR, "glide-array.sh")).read()
    assert "--array=" in text
    assert "SLURM_ARRAY_TASK_ID" in text


def test_longleaf_is_labeled_historical():
    text = open(os.path.join(SLURM_DIR, "longleaf-glide-gpu.sh")).read()
    assert "schrodinger/2021-4" in text
    assert "volta-gpu" in text
    assert "HISTORICAL" in text.upper() or "historical" in text


# ---- no-license guard ------------------------------------------------------

# Forbidden OSS-license identifiers, assembled from fragments so the tokens
# themselves never appear verbatim anywhere in the repo.
_IDS = ["M" + "IT", "Apa" + "che", "B" + "SD", "G" + "PL", "LG" + "PL", "M" + "PL"]
LICENSE_PATTERN = re.compile(r"\b(" + "|".join(_IDS) + r")\b")
GRANT_PATTERN = re.compile(r"permission is hereby " + "granted" + r"|licensed under", re.I)


def _repo_text_files():
    this_file = os.path.abspath(__file__)
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [
            d for d in dirnames
            if d not in {".git", "__pycache__", ".pytest_cache", ".egg-info"}
            and not d.endswith(".egg-info")
        ]
        for fn in filenames:
            full = os.path.join(dirpath, fn)
            if os.path.abspath(full) == this_file:
                continue  # this guard file names the patterns it forbids
            yield full


def test_no_license_file_anywhere():
    for path in _repo_text_files():
        base = os.path.basename(path).lower()
        assert base not in {"license", "license.txt", "license.md", "copying"}, path


def test_no_license_grant_in_sources():
    for path in _repo_text_files():
        if path.endswith((".png", ".jpg", ".zip", ".maegz", ".pdb")):
            continue
        try:
            text = open(path, encoding="utf-8").read()
        except (UnicodeDecodeError, OSError):
            continue
        assert not GRANT_PATTERN.search(text), f"license grant in {path}"
        assert not LICENSE_PATTERN.search(text), f"license identifier in {path}"

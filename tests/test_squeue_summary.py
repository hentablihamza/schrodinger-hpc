import os

import pytest

import squeue_summary as ss
from conftest import EXAMPLE_DIR

SAMPLE = os.path.join(EXAMPLE_DIR, "sacct_array_sample.txt")


def test_normalize_states():
    assert ss._normalize_state("COMPLETED") == "COMPLETED"
    assert ss._normalize_state("R") == "RUNNING"
    assert ss._normalize_state("CANCELLED by 42") == "FAILED"
    assert ss._normalize_state("TIMEOUT") == "FAILED"
    assert ss._normalize_state("WEIRD") == "OTHER"


def test_parse_skips_substeps_and_header():
    text = "JobID|State|Elapsed\n100_0|COMPLETED|1\n100_0.batch|COMPLETED|1\n100_1|FAILED|1\n"
    rows = ss.parse_states(text)
    assert rows == [("100_0", "COMPLETED"), ("100_1", "FAILED")]


@pytest.mark.skipif(not os.path.exists(SAMPLE), reason="sacct sample missing")
def test_summarize_sample():
    counts = ss.summarize(open(SAMPLE).read())
    # Sub-steps excluded; 10 array tasks counted once each.
    assert sum(counts.values()) == 10
    assert counts["COMPLETED"] == 4
    assert counts["FAILED"] == 2  # one FAILED + one TIMEOUT
    assert counts["RUNNING"] == 2
    assert counts["PENDING"] == 2


def test_squeue_whitespace_format():
    text = "JOBID ST\n55_0 R\n55_1 PD\n55_2 CD\n"
    counts = ss.summarize(text)
    assert counts["RUNNING"] == 1
    assert counts["PENDING"] == 1
    assert counts["COMPLETED"] == 1

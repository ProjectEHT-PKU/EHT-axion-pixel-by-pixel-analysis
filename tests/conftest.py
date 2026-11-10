"""Tests write stages from a working tree that may have uncommitted edits.

provenance.code_identity requires clean tracked files; for the whole test session it is
replaced by a fixed commit. tests/test_provenance.py checks the real function (captured
at import) with temporary git repositories.
"""
import pytest

from eht_axion import provenance


@pytest.fixture(scope="session", autouse=True)
def fixed_code_identity():
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(provenance, "code_identity", lambda root=None: {"git_commit": "0" * 40})
        yield

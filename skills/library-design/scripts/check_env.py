"""Environment guard for pegg library design. Run this BEFORE designing.

The pins this checks are load-bearing, not preferences. pegg's scoring models
(V3_model_full.pickle, V3_model_nopos.pickle, RF_PE.pkl) are pickles built
against scikit-learn 1.1.1; on a newer sklearn they raise an unpickling error
that does not name the pin, or silently mis-score. That takes out
ontarget_score and therefore PEGG2_Score and RF_Score -- discovered hours into
a run, not at import.

The pegg>=3.1.0 probe is lifted from ch_aml_library_design.ipynb cell 2:

    assert 'genomic_span' in inspect.signature(
        bystander.reverse_frame_anchor).parameters

Do not edit here without editing there -- these two must not drift. Under
3.0/3.0.1, reverse_frame_anchor had no genomic_span argument, and
silent_bystander=True on an INS/DEL/INDEL produced bystanders that were
non-synonymous in the real reading frame with no error raised.

Exits non-zero on anything that would silently produce a wrong library.
"""

import argparse
import inspect
import os
import sys


# (name, what it is for, required)
ENV_VARS = (
    ("LAB_DATA_ROOT", "reference genome and GENCODE .db", True),
    ("OPTIPRIME_SRC", "clone of alvin-hsu/optiprime-src (OptiPrime scoring)", False),
    ("OPTIPRIME_PYTHON", "the Python 3.11 interpreter OptiPrime runs in", False),
)


def _ok(msg):
    print("  ok    " + msg)


def _bad(msg):
    print("  FAIL  " + msg)


def check_python():
    v = sys.version_info
    got = "%d.%d.%d" % (v[0], v[1], v[2])
    if v[0] == 3 and v[1] in (9, 10):
        _ok("python " + got)
        return []
    _bad("python " + got + " -- pegg needs 3.9 or 3.10 "
         "(numpy<2 + scikit-learn==1.1.1 cap it there). "
         "conda env create -f envs/pegg_env.yml")
    return ["python"]


def check_pins():
    bad = []
    try:
        import sklearn
    except ImportError:
        _bad("scikit-learn is not installed -- needs exactly 1.1.1")
        bad.append("scikit-learn")
    else:
        if sklearn.__version__ == "1.1.1":
            _ok("scikit-learn " + sklearn.__version__)
        else:
            _bad("scikit-learn " + sklearn.__version__ + " -- needs exactly "
                 "1.1.1; the Azimuth and RF pickles will not load, so "
                 "PEGG2_Score and RF_Score will be missing or wrong")
            bad.append("scikit-learn")

    try:
        import numpy
    except ImportError:
        _bad("numpy is not installed -- needs >=1.21.5,<2")
        bad.append("numpy")
    else:
        major = int(numpy.__version__.split(".")[0])
        if major < 2:
            _ok("numpy " + numpy.__version__)
        else:
            _bad("numpy " + numpy.__version__ + " -- needs <2; 2.x breaks "
                 "binary compatibility with the sklearn 1.1.1 build and cyvcf2")
            bad.append("numpy")

    for mod in ("pandas", "Bio", "gffutils", "regex"):
        try:
            __import__(mod)
        except ImportError:
            _bad(mod + " is not importable (pegg needs it; gffutils is "
                 "required whenever silent bystanders are used)")
            bad.append(mod)
        else:
            _ok(mod + " importable")
    return bad


def check_pegg():
    bad = []
    try:
        import pegg
        from pegg import prime, bystander
    except Exception as exc:
        _bad("cannot import pegg: %s" % exc)
        return ["pegg"]

    ver = getattr(pegg, "__version__", "unknown")
    _ok("pegg " + str(ver) + " from " + os.path.dirname(pegg.__file__))

    if "silent_bystander" not in prime.run.__code__.co_varnames:
        _bad("this pegg has no silent_bystander argument -- it is older than "
             "3.0. pip install --upgrade 'pegg>=3.1.0'")
        bad.append("pegg")
        return bad

    params = inspect.signature(bystander.reverse_frame_anchor).parameters
    if "genomic_span" in params:
        _ok("pegg >= 3.1.0 (reverse_frame_anchor takes genomic_span) -- "
            "silent bystanders on indels are frame-correct")
    else:
        _bad("pegg is older than 3.1.0: reverse_frame_anchor has no "
             "genomic_span argument, so silent_bystander=True on INS/DEL/"
             "INDEL produces NON-SYNONYMOUS bystanders with no error. "
             "pip install --upgrade 'pegg>=3.1.0' and REGENERATE any indel "
             "library designed with this install.")
        bad.append("pegg")
    return bad


def check_env_vars(need_optiprime):
    bad = []
    for name, what, required in ENV_VARS:
        val = os.environ.get(name)
        required = required or (need_optiprime and name.startswith("OPTIPRIME"))
        if not val:
            if required:
                _bad(name + " is not set (" + what + ") -- see docs/data-layout.md")
                bad.append(name)
            else:
                print("  --    " + name + " unset (" + what + "); only needed "
                      "for OptiPrime scoring")
            continue
        exists = os.path.isdir(val) or os.path.isfile(val)
        if exists:
            _ok(name + " = " + val)
        else:
            _bad(name + " = " + val + " -- path does not exist")
            bad.append(name)
    return bad


def check_optiprime():
    try:
        from pegg import optiprime as _op
    except Exception as exc:
        _bad("cannot import pegg.optiprime: %s" % exc)
        return ["optiprime"]
    if _op.optiprime_available():
        _ok("OptiPrime configured and reachable (default group %s)"
            % _op.DEFAULT_GROUP)
        return []
    _bad("OptiPrime is not reachable. Note a Jupyter kernel does NOT inherit "
         "the shell environment -- set OPTIPRIME_SRC / OPTIPRIME_PYTHON with "
         "os.environ[...] inside the notebook, or pass optiprime_src= / "
         "optiprime_python= to prime.run().")
    return ["optiprime"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--optiprime", action="store_true",
                    help="also require a working OptiPrime install")
    args = ap.parse_args()

    print("python / pins")
    bad = check_python()
    bad += check_pins()
    print("pegg")
    bad += check_pegg()
    print("environment")
    bad += check_env_vars(args.optiprime)
    if args.optiprime:
        bad += check_optiprime()

    print("")
    if bad:
        raise SystemExit("NOT SAFE TO DESIGN -- fix: " + ", ".join(sorted(set(bad))))
    print("environment ok")


if __name__ == "__main__":
    main()

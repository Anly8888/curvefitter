#!/usr/bin/env python3
"""Regenerate oracle fixtures using pinned, unmodified ImageJ Java sources (JDK 21)."""
import argparse
import hashlib
from pathlib import Path
import subprocess
import tempfile

COMMIT = "888e13d30cf0fa7c13f601fe929fdb3f940cafea"
HASHES = {
    "ij/measure/CurveFitter.java": "b71fc14eb53560bba71e8921653b0dcbafaf9f06d43edf4e7cd35f571365aced",
    "ij/measure/Minimizer.java": "9b850720f808716155efdd9f57eaff139e8960276c71c7e3b87dea254a5dd776",
    "ij/util/IJMath.java": "18dd0663e39ecb5c354886cd2bce1ac60489bd390d24bcaee8ad1c4cc0a55207",
}


def run(*args, **kwargs):
    return subprocess.run(args, check=True, **kwargs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, help="Existing ImageJ checkout at the pinned commit")
    parser.add_argument("--check", action="store_true", help="Compare generated output without writing files")
    args = parser.parse_args()
    tests = Path(__file__).resolve().parent
    with tempfile.TemporaryDirectory(prefix="curvefitter-oracle-") as temp:
        temp = Path(temp)
        source = args.source.resolve() if args.source else temp / "imagej"
        if not args.source:
            run("git", "init", "-q", str(source))
            run("git", "-C", str(source), "remote", "add", "origin", "https://github.com/imagej/ImageJ.git")
            run("git", "-C", str(source), "fetch", "--depth=1", "origin", COMMIT)
            run("git", "-C", str(source), "checkout", "--detach", "FETCH_HEAD")
        head = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
        if head != COMMIT:
            raise RuntimeError(f"Expected ImageJ commit {COMMIT}, found {head}")
        for path, expected in HASHES.items():
            actual = hashlib.sha256((source / path).read_bytes()).hexdigest()
            if actual != expected:
                raise RuntimeError(f"Unmodified upstream source required: {path}")
        classes = temp / "classes"
        classes.mkdir()
        compiler = ("java", "-m", "jdk.compiler/com.sun.tools.javac.Main")
        run(*compiler, "-nowarn", "-d", str(classes), "ij/ImageJ.java", cwd=source)
        run(*compiler, "-cp", str(classes), "-d", str(classes), str(tests / "javarunner.java"))
        output = subprocess.check_output(["java", "-Djava.awt.headless=true", "-cp", str(classes), "JavaRunner"])
        target = tests / "java_reference.jsonl"
        if args.check:
            if output != target.read_bytes():
                raise RuntimeError("Java reference fixtures differ; inspect before updating")
        else:
            target.write_bytes(output)
        print(f"Verified {len(output.splitlines())} Java reference cases at {COMMIT}")


if __name__ == "__main__":
    main()

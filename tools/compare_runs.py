"""Compare CLI results and alternating timings of two builds (no dependencies)."""
import argparse
import json
from pathlib import Path
import random
import re
import statistics
import subprocess
import tempfile
import time


def run(binary, args, timeout=180):
    start = time.perf_counter()
    result = subprocess.run([str(binary), *args], capture_output=True, text=True,
                            timeout=timeout)
    return result, time.perf_counter() - start


def stable_output(text):
    # Only benchmark timing fields vary; polynomial and summary lines stay exact.
    return re.sub(r"(Time: |time: )\d+(?:\.\d+)?(?:e[+-]?\d+)?", r"\1<TIME>", text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("repeats must be positive")
    baseline, candidate = args.baseline.resolve(), args.candidate.resolve()
    rng = random.Random(20261011)
    lines = []
    for n in range(2, 10):
        for _ in range(40):
            x, w = list(range(n)), list(range(n))
            rng.shuffle(x)
            rng.shuffle(w)
            lines.append("".join(map(str, x)) + " " + "".join(map(str, w)))
        for _ in range(10):
            w = list(range(n))
            rng.shuffle(w)
            lines.append("".join(map(str, range(n))) + " " + "".join(map(str, w)))
    with tempfile.TemporaryDirectory(dir=candidate.parent) as tmp:
        pairs = Path(tmp) / "pairs.txt"
        pairs.write_text("\n".join(lines), encoding="utf-8")
        commands = [["--verify", "verify.txt"], ["--file", str(pairs)],
                    ["0123", "2301"], ["--simple", "2301"],
                    ["--random", "8", "20", "42"], ["--cyclic", "10"]]
        for command in commands:
            old, _ = run(baseline, command)
            new, _ = run(candidate, command)
            if (old.returncode, stable_output(old.stdout), old.stderr) != (new.returncode, stable_output(new.stdout), new.stderr):
                raise AssertionError(f"Output changed: {command}\nOLD: {old}\nNEW: {new}")
        print(f"Exact CLI output comparison passed ({len(lines)} generated pairs + verify.txt).", flush=True)
    benchmarks = [["--speedtest2", "10", "500", "42"],
                  ["--speedtest2", "11", "100", "42"],
                  ["0123456789a", "789a3456012"],
                  ["--verify", "verify.txt"]]
    report = []
    for command in benchmarks:
        measurements = [[], []]
        outputs = []
        for binary in (baseline, candidate):
            result, _ = run(binary, command)  # Fresh calculator; warm filesystem/code pages.
            if result.returncode:
                raise AssertionError(result.stderr)
            outputs.append(stable_output(result.stdout))
        if outputs[0] != outputs[1]:
            raise AssertionError(f"Benchmark output changed: {command}")
        for i in range(args.repeats):
            for index in ([0, 1] if i % 2 == 0 else [1, 0]):
                result, elapsed = run((baseline, candidate)[index], command)
                if result.returncode:
                    raise AssertionError(result.stderr)
                if stable_output(result.stdout) != outputs[index]:
                    raise AssertionError(f"Inconsistent benchmark output: {command}")
                measurements[index].append(elapsed)
        old, new = map(statistics.median, measurements)
        item = dict(command=command, baseline_seconds=old, candidate_seconds=new,
                    ratio=new / old, samples=measurements)
        report.append(item)
        print(json.dumps(item), flush=True)
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

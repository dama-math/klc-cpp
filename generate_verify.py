"""Generate regression fixtures. Expected values come from the chosen binary,
so use tests/core_tests.cpp for independent mathematical validation.
"""
import argparse
from pathlib import Path
import random
import subprocess

def to_one_line(perm):
    res = ""
    for v in perm:
        if v < 10:
            res += str(v)
        else:
            res += chr(ord('a') + v - 10)
    return res

def run_klc(binary, x_str, w_str, timeout):
    try:
        result = subprocess.run([str(binary), x_str, w_str], capture_output=True,
                                text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None
    if result.returncode != 0:
        raise RuntimeError(f"Calculator failed: {result.stderr.strip()}")
    parts = result.stdout.split()
    if len(parts) != 3 or parts[:2] != [x_str, w_str]:
        raise RuntimeError(f"Unexpected calculator output: {result.stdout!r}")
    tuple(map(int, parts[2].split(",")))
    return parts[2]


def is_3412_avoiding(w):
    n = len(w)
    for i in range(n):
        for j in range(i + 1, n):
            if w[i] < w[j]:
                for k in range(j + 1, n):
                    for l in range(k + 1, n):
                        if w[k] < w[l] < w[i]:
                            return False
    return True

# (size range, quotas, use identity x, require 3412 avoidance, minimum P(1), title)
GROUPS = [
    (5, 8, {"zero": 2, "trivial": 5, "nontrivial": 8}, False, False, 0,
     "Test cases for S_5 to S_8 (15 cases: 8 non-trivial, 5 trivial, 2 zero)"),
    (9, 11, {"zero": 2, "trivial": 2, "nontrivial": 8}, False, False, 0,
     "Test cases for S_9 to S_11 (12 cases: 8 non-trivial, 2 trivial, 2 zero)"),
    (11, 11, {"large": 3}, False, False, 50,
     "Test cases for n = 11 with medium complexity (3 cases, P(1) >= 50)"),
    (12, 12, {"large": 5}, True, False, 1001,
     "Test cases for n = 12 with large KL polynomials (5 cases, P(1) > 1000)"),
    (12, 12, {"avoiding": 5}, True, True, 0,
     "Test cases for n = 12 with 3412-avoiding permutations w (5 cases, Lascoux optimization)"),
]


def generate_group(binary, rng, group, max_attempts, timeout):
    minimum, maximum, quota, identity_x, avoiding, threshold, title = group
    remaining = quota.copy()
    lines = ["# " + title]
    recorded = set()
    for _ in range(max_attempts):
        n = rng.randint(minimum, maximum)
        w, x = list(range(n)), list(range(n))
        rng.shuffle(w)
        if avoiding and not is_3412_avoiding(w):
            continue
        if not identity_x:
            rng.shuffle(x)
        pair = to_one_line(x), to_one_line(w)
        if pair in recorded:
            continue
        poly = run_klc(binary, *pair, timeout)
        if poly is None:
            continue
        if avoiding:
            category = "avoiding" if poly != "0" else "zero"
        elif threshold:
            category = "large" if sum(map(int, poly.split(","))) >= threshold else "other"
        else:
            category = "zero" if poly == "0" else "nontrivial" if "," in poly else "trivial"
        if remaining.get(category, 0):
            lines.append(f"{pair[0]} {pair[1]} {poly}")
            recorded.add(pair)
            remaining[category] -= 1
        if not any(remaining.values()):
            return lines
    raise RuntimeError(f"Generation attempt limit reached for {title}; remaining: {remaining}")


def generate(binary, output, seed, timeout, max_attempts, groups=GROUPS):
    rng = random.Random(seed)
    lines = []
    for group in groups:
        if lines:
            lines.append("")
        lines.extend(generate_group(binary, rng, group, max_attempts, timeout))
    # Failed generation never replaces a previous fixture.
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    count = sum(bool(line) and not line.startswith("#") for line in lines)
    print(f"Generated {output} with {count} test cases across {len(groups)} groups.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=Path("klc_cpp.exe") if Path("klc_cpp.exe").exists() else Path("klc_cpp"))
    parser.add_argument("--output", type=Path, default=Path("verify.generated.txt"))
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--timeout", type=float, default=30)
    parser.add_argument("--max-attempts", type=int, default=100000)
    parser.add_argument("--group", type=int, choices=range(1, 6))
    args = parser.parse_args()
    if args.timeout <= 0 or args.max_attempts <= 0:
        parser.error("timeout and max-attempts must be positive")
    groups = GROUPS if args.group is None else [GROUPS[args.group - 1]]
    generate(args.binary.resolve(strict=True), args.output, args.seed, args.timeout, args.max_attempts, groups)


if __name__ == "__main__":
    main()

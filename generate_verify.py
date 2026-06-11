import random
import subprocess
import sys

def to_one_line(perm):
    res = ""
    for v in perm:
        if v < 10:
            res += str(v)
        else:
            res += chr(ord('a') + v - 10)
    return res

def run_klc(x_str, w_str):
    result = subprocess.run(["./klc_cpp", x_str, w_str], capture_output=True, text=True)
    if result.returncode != 0:
        return None
    out = result.stdout.strip()
    if not out:
        return None
    lines = out.split('\n')
    parts = lines[-1].split()
    if len(parts) < 3:
        return None
    return parts[2]

verify_lines = []

# --- 5 to 8 group (15 cases: 8 non-trivial, 5 trivial, 2 zero) ---
zero_needed = 2
trivial_needed = 5
nontrivial_needed = 8

while zero_needed > 0 or trivial_needed > 0 or nontrivial_needed > 0:
    n = random.randint(5, 8)
    w = list(range(n))
    random.shuffle(w)
    x = list(range(n))
    random.shuffle(x)
    
    w_str = to_one_line(w)
    x_str = to_one_line(x)
    
    poly = run_klc(x_str, w_str)
    if not poly:
        continue
        
    if poly == "0":
        if zero_needed > 0:
            verify_lines.append(f"{x_str} {w_str} {poly}")
            zero_needed -= 1
    else:
        is_trivial = ("," not in poly)
        if is_trivial and trivial_needed > 0:
            verify_lines.append(f"{x_str} {w_str} {poly}")
            trivial_needed -= 1
        elif not is_trivial and nontrivial_needed > 0:
            verify_lines.append(f"{x_str} {w_str} {poly}")
            nontrivial_needed -= 1

# --- 9 to 11 group (12 cases: 8 non-trivial, 2 trivial, 2 zero) ---
zero_needed = 2
trivial_needed = 2
nontrivial_needed = 8

while zero_needed > 0 or trivial_needed > 0 or nontrivial_needed > 0:
    n = random.randint(9, 11)
    w = list(range(n))
    random.shuffle(w)
    x = list(range(n))
    random.shuffle(x)
    
    w_str = to_one_line(w)
    x_str = to_one_line(x)
    
    poly = run_klc(x_str, w_str)
    if not poly:
        continue
        
    if poly == "0":
        if zero_needed > 0:
            verify_lines.append(f"{x_str} {w_str} {poly}")
            zero_needed -= 1
    else:
        is_trivial = ("," not in poly)
        if is_trivial and trivial_needed > 0:
            verify_lines.append(f"{x_str} {w_str} {poly}")
            trivial_needed -= 1
        elif not is_trivial and nontrivial_needed > 0:
            verify_lines.append(f"{x_str} {w_str} {poly}")
            nontrivial_needed -= 1

# --- n = 11 group (3 cases, sum >= 50) ---
needed = 3
while needed > 0:
    n = 11
    w = list(range(n))
    random.shuffle(w)
    x = list(range(n))
    random.shuffle(x)
    
    w_str = to_one_line(w)
    x_str = to_one_line(x)
    
    poly = run_klc(x_str, w_str)
    if not poly or poly == "0":
        continue
        
    coeffs = [int(c) for c in poly.split(",")]
    if sum(coeffs) >= 50:
        verify_lines.append(f"{x_str} {w_str} {poly}")
        needed -= 1

with open("verify.txt", "w") as f:
    for line in verify_lines:
        f.write(line + "\n")

print(f"Generated verify.txt with {len(verify_lines)} cases.")

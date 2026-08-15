import os
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

klc_binary = "./klc_cpp.exe" if os.path.exists("./klc_cpp.exe") else "./klc_cpp"

def run_klc(x_str, w_str):
    result = subprocess.run([klc_binary, x_str, w_str], capture_output=True, text=True)
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

verify_lines = []

# ==============================================================================
# Group 1: S_5 to S_8 (15 cases: 8 non-trivial, 5 trivial, 2 zero)
# ==============================================================================
verify_lines.append("# Test cases for S_5 to S_8 (15 cases: 8 non-trivial, 5 trivial, 2 zero)")
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

# ==============================================================================
# Group 2: S_9 to S_11 (12 cases: 8 non-trivial, 2 trivial, 2 zero)
# ==============================================================================
verify_lines.append("\n# Test cases for S_9 to S_11 (12 cases: 8 non-trivial, 2 trivial, 2 zero)")
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

# ==============================================================================
# Group 3: n = 11 (3 cases, sum of coefficients >= 50)
# ==============================================================================
verify_lines.append("\n# Test cases for n = 11 with medium complexity (3 cases, P(1) >= 50)")
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

# ==============================================================================
# Group 4: n = 12 (5 cases, P(1) > 1000)
# ==============================================================================
verify_lines.append("\n# Test cases for n = 12 with large KL polynomials (5 cases, P(1) > 1000)")
needed = 5
while needed > 0:
    n = 12
    w = list(range(n))
    random.shuffle(w)
    x = list(range(n))
    
    w_str = to_one_line(w)
    x_str = to_one_line(x)
    
    poly = run_klc(x_str, w_str)
    if not poly or poly == "0":
        continue
        
    coeffs = [int(c) for c in poly.split(",")]
    if sum(coeffs) > 1000:
        verify_lines.append(f"{x_str} {w_str} {poly}")
        needed -= 1

# ==============================================================================
# Group 5: n = 12 with 3412-avoiding permutation w (5 cases, testing Lascoux's method)
# ==============================================================================
verify_lines.append("\n# Test cases for n = 12 with 3412-avoiding permutations w (5 cases, Lascoux optimization)")
needed = 5
while needed > 0:
    n = 12
    w = list(range(n))
    random.shuffle(w)
    if not is_3412_avoiding(w):
        continue
        
    x = list(range(n))
    w_str = to_one_line(w)
    x_str = to_one_line(x)
    
    poly = run_klc(x_str, w_str)
    if not poly or poly == "0":
        continue
        
    verify_lines.append(f"{x_str} {w_str} {poly}")
    needed -= 1

with open("verify.txt", "w") as f:
    for line in verify_lines:
        f.write(line + "\n")

num_cases = sum(1 for line in verify_lines if line and not line.startswith("#") and not line.startswith("\n"))
print(f"Generated verify.txt with {num_cases} test cases across 5 groups.")

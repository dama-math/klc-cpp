import random
import subprocess
import string
import sys

def to_one_line(perm):
    res = ""
    for v in perm:
        if v < 10:
            res += str(v)
        else:
            res += chr(ord('a') + v - 10)
    return res

with open("verify.txt", "w") as f:
    for _ in range(10):
        n = random.randint(5, 11)
        w = list(range(n))
        random.shuffle(w)
        x = list(range(n))
        
        w_str = to_one_line(w)
        x_str = to_one_line(x)
        
        result = subprocess.run(["./klc_cpp", x_str, w_str], capture_output=True, text=True)
        if result.returncode == 0:
            out = result.stdout.strip()
            # It should output exactly 'x_str w_str poly' without any internal newlines
            if '\n' not in out and out:
                f.write(out + '\n')
            else:
                # If there's some extra output, just take the last line assuming it's the result
                lines = out.split('\n')
                f.write(lines[-1] + '\n')

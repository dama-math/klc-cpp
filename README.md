# Type-A Kazhdan-Lusztig Polynomial Calculator (`klc-cpp`)

A C++ program for calculating Kazhdan-Lusztig (KL) polynomials $P_{x,w}(q)$ for permutations in symmetric groups $S_n$.

Based on Gregory S. Warrington's original C implementation ([Warrington's KLC page](https://gswarrin.w3.uvm.edu/research/klc/klc.html)).

---

## Features

- **Recurrence solver:** C++ port of Warrington's `klc` code.
- **Vexillary permutations:** Uses Lascoux's algorithm for 3412-avoiding (vexillary) permutations to speed up applicable cases.
- **Compile-time integer packing:**
  - **64-bit (default):** Packs elements into a `uint64_t` (4 bits/element), supporting $n \le 16$.
  - **128-bit:** Packs elements into `__uint128_t` (5 bits/element) when compiled with `BIG_N=1`, supporting $n \le 24$.
- **AVX2 inversion counting:** Uses SIMD vectorization for Coxeter length calculation when the host CPU supports it.
- **Caching:** Previously computed polynomials and $\mu$ coefficients are cached to speed up recursive lookups.

---

## Permutation Notation & Math Constraints

### Permutation representation

Permutations are in **0-indexed one-line notation** over $\{0, 1, \dots, n-1\}$.

Every element from $0$ to $n-1$ must appear exactly once.

- Valid ($S_4$): `0213`
- Invalid ($S_4$): `1234` (element `4` is out of bounds; `0` is missing)

For elements 10 and above, use lowercase letters: `a` = 10, `b` = 11, ..., `n` = 23.

- Example ($S_{12}$ identity): `0123456789ab`

### Mathematical constraints

- $P_{x,w}(q)$ is computed only when $x \le w$ in the Bruhat order. If $x \not\le w$, the program returns `0`.
- Output is a comma-separated list of coefficients starting from degree 0. For example, `1,3,2` represents $1 + 3q + 2q^2$.

---

## Build Instructions

Requires a C++20 compliant compiler (`clang++` or `g++`) and `make`.

| Target | Command | Max $n$ |
| :--- | :--- | :---: |
| Standard | `make` | $\le 16$ |
| Large-N | `make BIG_N=1` | $\le 24$ |

To use `g++` instead of the default `clang++`:

```bash
make COMPILER=gcc
```

To clean build artifacts:

```bash
make clean
```

---

## Usage

```bash
./klc_cpp [options]
./klc_cpp x_permutation w_permutation
```

| Option | Argument | Description |
| :--- | :--- | :--- |
| (none) | `x_perm w_perm` | Calculates $P_{x,w}(q)$. |
| `--help` | — | Displays help. |
| `--simple` | `w_perm` | Calculates $P_{e,w}(1)$, where $e$ is the identity. |
| `--speedtest` | — | Runs a benchmark (all KL polynomials for $n \le 9$). |
| `--random` | `n [num] [seed]` | Runs `num` random KL calculations in $S_n$. |
| `--file` | `filename` | Processes permutation pairs from a file. |
| `--verify` | `filename` | Verifies calculations against expected results in a file. |

### Examples

Calculate $P_{x,w}(q)$ for a specific pair in $S_4$:

```bash
./klc_cpp 0123 2301
# Output: 0123 2301 1,1
# (Representing 1 + q)
```

Calculate $P_{e,w}(1)$:

```bash
./klc_cpp --simple 2301
# Output: 2
```

Verify against a file:

```bash
./klc_cpp --verify verify.txt
```

**File formats:**

- `--file`: Each line contains a space-separated pair `x_permutation w_permutation`.
- `--verify`: Each line contains `x_permutation w_permutation expected_coefficients` (e.g. `0213 2301 1,1`).

---

## References

- Gregory S. Warrington, original C implementation ([official webpage](https://gswarrin.w3.uvm.edu/research/klc/klc.html))
- A. Lascoux, *Polynômes de Kazhdan-Lusztig pour les variétés de Schubert vexillaires*, C. R. Acad. Sci. Paris Sér. I Math. 321 (1995), 667–670.
- Sara Billey and V. Lakshmibai, *Singular Loci of Schubert Varieties*, Vol. 182. Springer Science & Business Media, 2000 (Section 6.3.28).

---

## License

GNU General Public License v3.0 — see [LICENSE](LICENSE).
# Kazhdan-Lusztig Polynomial Calculator (klc-cpp-release)

This is a high-performance C++ program for calculating Kazhdan-Lusztig (KL) polynomials for permutations in symmetric groups.

This calculator is based on the C implementation by Gregory S. Warrington (available at https://gswarrin.w3.uvm.edu/research/klc/klc.html).

The calculator utilizes various optimizations to speed up the computation, such as:
- Lascoux's algorithm for 3412-free (vexillary) permutations.
  - Original paper: A. Lascoux, *Polynômes de Kazhdan-Lusztig pour les variétés de Schubert vexillaires*, C. R. Acad. Sci. Paris Sér. I Math. 321 (1995), 667–670.
  - See also: Sara Billey and V. Lakshmibai, *Singular Loci of Schubert Varieties* (Section 6.3.28).
- Efficient memory layouts and bit manipulation.

## Build Instructions

This project features a **compile-time hybrid architecture** to provide maximum performance for small permutations and safe computation for large ones.

**For maximum performance ($n \le 16$):**
```bash
make
```
This builds the default 64-bit mode, which is extremely fast and heavily optimized for cache locality.

**For large permutations ($n \le 24$):**
```bash
make BIG_N=1
```
This builds the 128-bit mode, which safely packs elements into native `__uint128_t` integers to support huge permutation spaces without memory overflow.

## Usage

The program provides a command-line interface with several modes of operation.

```bash
Usage: ./klc_cpp [options]
Options:
  --speedtest                       Speed test KL polynomial calculations
  --random n [num] [seed]           Perform random KL calculations for S_n
  --verify filename                 Verify KL polynomial calculations against a file
  --cyclic n                        Run a fast cyclic permutation test for S_n
  --file filename                   Process permutations from a file
  --simple w_permutation            Calculate KL polynomial P_e_w(1) for given permutation w
  x_permutation w_permutation       Calculate KL polynomial for given permutations
```

### Examples

**Calculate the KL polynomial for a specific pair of permutations ($x$ and $w$):**
```bash
./klc_cpp 0123 3210
```

**Calculate $P_{e,w}(1)$ for a specific permutation $w$:**
```bash
./klc_cpp --simple 3210
```

**Perform a speed test (calculates KL polynomials for $n \le 9$):**
```bash
./klc_cpp --speedtest
```

**Perform random KL calculations in $S_{10}$:**
```bash
./klc_cpp --random 10 5
```

**Process a file containing pairs of permutations:**
```bash
./klc_cpp --file input.txt
```
*(The input file should contain one pair of space-separated permutations per line.)*

**Verify polynomial calculations against a file:**
```bash
./klc_cpp --verify verify.txt
```
*(The verification file should contain space-separated `x w expected_polynomial` per line.)*

## Permutation Notation

Permutations are represented in one-line notation. For elements $10$ and above, alphabetical characters are used (`a` for $10$, `b` for $11$, etc.). For example, the identity permutation in $S_{12}$ is represented as `0123456789ab`.

> [!IMPORTANT]
> **Maximum Length:** The permutation length $n$ must be **$\le 16$** by default, or **$\le 24$** when compiled with `make BIG_N=1`. 
> - In **64-bit mode** (default), elements are packed into a `uint64_t` using 4 bits per element, supporting up to $n=16$ with extreme cache-efficiency.
> - In **128-bit mode** (`BIG_N=1`), elements are packed into a `__uint128_t` using 5 bits per element, safely supporting up to $n=24$.

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

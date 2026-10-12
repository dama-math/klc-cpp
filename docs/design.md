# Structure and invariants

`main.cpp` only enters the CLI. `src/Cli.cpp` parses complete integer arguments,
validates argument counts, and translates errors into process exit status.
`src/Commands.cpp` handles file I/O, formatting, random inputs and benchmark clocks.
The calculator and Lascoux implementation do not print library diagnostics.

`Permutation` is a fixed-size value type. All public constructors validate size,
range and uniqueness before the data can reach bit operations or inverse indexing.
`identity()` and simple transpositions keep the inversion cache consistent.
`at()` checks an external index; `operator[]` has a valid-index precondition and a
debug assertion, like other indexed accesses in the internal hot loops.
Mutable storage and cache/size updates are restricted to internal operations.
`detail/PermutationAccess.hpp` is reserved for reduction and candidate enumeration;
its commit operation checks both permutation validity and cached length in debug
builds. These internal operations must restore invariants before returning.

`Polynomial` uses a fixed array (64 coefficients in standard mode, 139 in
large mode) and degree -1 for zero. Shifted accumulation is shared by addition and
subtraction, supports self-aliasing, checks capacity and performs checked signed
arithmetic. Multiplication and evaluation also detect overflow. Polynomial
arithmetic is in-place: if arithmetic throws, callers should discard the partial
result. Constructors and multiplication do not commit a partial result.

`KLPolynomialCalculator` owns its two memoization tables. Pair reduction and
inverse canonicalization occur before lookup. KL and mu recurrences deliberately
remain separate because their correction order and early-exit rules differ.
The same canonicalization helper defines their key ordering.

`detail/IntervalSearch.hpp` contains one interval enumeration algorithm. Its rank
backend is selected once per search: SSE4.1 for sizes <=16, scalar otherwise.
`detail/RankOperations.hpp` contains the only SIMD intrinsics, aligned lookup
tables and lane-width assumptions. Rank counts and differences fit in int8_t
because the supported maximum size is 24. Empty and singleton cases return
before computing any descent masks or shifts. Both backends use identical
candidate order and callback termination semantics.

`detail/FlatMap.hpp` is an open-addressed cache specialized in use to packed keys
and trivially copied results. Rehashing is prepared in a temporary table so an
allocation failure preserves the original table. A value that aliases an entry
is copied before growth. Lookup pointers must not be held across an insertion
that can grow the cache. Clear retains storage for reuse.

Tests use independent prefix-rank and unoptimized recurrence oracles, and compare
enumeration against all candidates in S_5. `verify.txt` is a frozen
regression fixture. `tools/compare_runs.py` compares stdout, stderr and exit status
for valid CLI commands, then alternates fresh processes and reports median wall
times. Timing fields are excluded from equality checks; the polynomial and other
text are compared exactly. Hardware, compiler flags and sample arrays should be
recorded alongside any performance conclusion.

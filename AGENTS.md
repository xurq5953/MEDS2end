# AGENTS.md

## Project purpose

This repository starts from the MEDS reference implementation and is being developed into a modified signature scheme.

The current development direction introduces:

- three equal dimensions;
- trilinear-form slice representations;
- corank-one points;
- canonical-form computation (`CF`);
- `BuildUVW`;
- `DiagonalNormalize`;
- new finite-field and square-matrix primitives required by those algorithms.

This is not a generic cleanup project. Changes must preserve the exact mathematical meaning of the cryptographic construction.

The repository is currently in a transition state:

- existing KeyGen, Sign, Verify, serialization, and the legacy three-matrix action may still represent the pre-CF MEDS2endGen design;
- the new low-level field and linear-algebra operators are intended to support the corank/CF design;
- do not silently rewrite upper-layer behavior while working on low-level primitives;
- upper-layer migration must be performed through an explicit plan and independently reviewed commits.

Preserve working MEDS infrastructure where it remains compatible, but do not retain legacy behavior when it conflicts with the latest authoritative algorithm specification.

---

## Instruction and specification priority

Use the following priority order:

1. The user’s current task prompt.
2. The latest user-provided algorithm specification introducing Corank and CF.
3. Any current repository specification file explicitly identified by the user.
4. Current implementation plans supplied for this repository.
5. This `AGENTS.md`.
6. The original MEDS documentation, used only as background for legacy code structure, serialization, seed handling, and notation.
7. Existing source-code conventions.

When two sources disagree, follow the higher-priority source.

Do not assume that an older `meds2end.tex`, an original MEDS paper, or an existing C function still describes the current target algorithm.

Before modifying any of the following, read the latest algorithm specification and inspect all affected call sites:

- KeyGen;
- Sign;
- Verify;
- trilinear-form actions;
- `BuildUVW`;
- `DiagonalNormalize`;
- `CF`;
- public-key format;
- secret-key format;
- signature format;
- challenge generation;
- response generation;
- serialization and hashing.

If the latest algorithm document is not present in the repository, do not reconstruct missing details from memory. Use the version supplied with the task or ask for the exact specification.

---

## Dimension and notation rules

The new algorithm uses equal dimensions:

```text
m = n = l
```

The current C code may still use historical MEDS macro names:

```text
MEDS_n
MEDS_n
MEDS_n
```

For the current equal-dimension parameter set, these are expected to satisfy:

```text
MEDS_n == MEDS_n
MEDS_n == MEDS_n
```

Use macros instead of hard-coded numeric dimensions.

For new corank/CF linear-algebra code:

- treat all principal matrices as `n x n` square matrices;
- treat vectors as length-`n` column vectors;
- do not introduce a general rectangular-matrix abstraction unless explicitly required;
- the existing rectangular helper may remain for legacy MEDS compatibility.

The current parameter set is expected to use dimension 14, but implementation code must not hard-code `14`.

---

## Repository structure

Inspect the actual repository before editing. Expected core files include:

- `meds.c`: implementation of `crypto_sign_keypair`, `crypto_sign`, and `crypto_sign_open`;
- `meds.h`: scheme-facing declarations;
- `params.h`: field type, parameters, dimensions, and byte-size macros;
- `api.h`: NIST/PQC API sizes and algorithm name;
- `field.c`, `field.h`: finite-field scalar arithmetic;
- `matrixmod.c`, `matrixmod.h`: non-elimination matrix and vector operations plus legacy systematic-form helpers;
- `matrixelim.c`, `matrixelim.h`: pivot-aware variable-time elimination, rank, inverse, kernels, and span tests;
- `triform.c`, `triform.h`: trilinear-form slice access, derived matrices, evaluation, corank-one wrappers, and the independent pullback action;
- `util.c`, `util.h`: XOF, random field elements, random matrices, challenge parsing, and legacy transformation helpers;
- `bitstream.c`, `bitstream.h`: field-element serialization;
- `fips202.c`, `fips202.h`: SHAKE/XOF implementation;
- `randombytes.c`, `randombytes.h`: randomness interface;
- `bench.c`: benchmark target;
- `PQCgenKAT_sign.c`: KAT generator;
- `Makefile`;
- `CMakeLists.txt`.

Detailed functional and stress tests should normally live only on a dedicated test branch.

---

# Finite-field layer

## Ownership

Finite-field scalar arithmetic belongs in:

```text
field.h
field.c
```

Do not move field arithmetic back into `matrixmod.c`.

`field.h` currently provides or is expected to provide:

```c
GF_add()
GF_sub()
GF_neg()
GF_mul()
GF_inv()
GF_inv_checked()
GF_batch_inv()
```

## Canonical representation

Every ordinary field element must be represented canonically:

```text
0 <= a < MEDS_p
```

All public field helpers assume canonical inputs unless their interface explicitly states otherwise.

Do not pass unchecked serialized values directly into arithmetic. Decoders must reject values greater than or equal to `MEDS_p`.

## Arithmetic rules

Use the existing field helpers instead of manually duplicating modular arithmetic.

Preferred usage:

```c
Fq c = GF_add(a, b);
Fq d = GF_sub(a, b);
Fq e = GF_mul(a, b);
```

Avoid open-coded expressions such as:

```c
(a + b) % MEDS_p
(a - b + MEDS_p) % MEDS_p
```

inside new matrix algorithms unless a measured optimization has been independently verified.

## Inversion rules

`GF_inv()` has a nonzero-input precondition:

```text
a != 0
```

Do not call `GF_inv(0)`.

When zero is a possible input, use:

```c
GF_inv_checked(&out, a)
```

and propagate failure.

This is mandatory in normalization code where a zero anchor must make the operation fail.

## Batch inversion

`GF_batch_inv()` is intended for collections of known nonzero field elements, such as normalization anchors.

Expected behavior:

- `count == 0` succeeds without accessing elements;
- any zero input causes failure;
- all nonzero inputs are inverted with one field inversion and linear additional work.

Current aliasing rule:

```text
out and in must not overlap
```

Do not call:

```c
GF_batch_inv(values, values, count);
```

unless the implementation is explicitly changed and tested for in-place operation.

---

# Matrix representation

## Element type

The only finite-field element type in the project is:

```c
Fq
```

Matrices, vectors, kernels, bases, coefficients, and trilinear-form slices are represented as contiguous arrays of `Fq`.

Do not introduce another scalar, vector, matrix, or finite-field element alias without an explicit migration plan.

Do not reintroduce finite-field element aliases.

## Storage order

Matrices use row-major storage.

For an `n x n` matrix:

```c
A[row * n + col]
```

For legacy rectangular helpers:

```c
A[row * number_of_columns + col]
```

Do not transpose the storage convention.

## Vector convention

Vectors are mathematical column vectors stored as contiguous arrays:

```c
Fq v[n];
```

A matrix-vector product computes:

```text
out = A v
```

The transpose-vector helper computes:

```text
out = A^T v
```

## Basis-vector layout

For span and independence helpers, a list of vectors is stored as consecutive vectors:

```text
basis[0 * n ... 1 * n - 1]   = vector 0
basis[1 * n ... 2 * n - 1]   = vector 1
...
```

When converted to a matrix, these vectors are inserted as columns.

Do not reinterpret `basis` as row vectors.

---

# Non-elimination matrix layer

## Ownership

Operations that do not require pivot search or Gaussian elimination belong in:

```text
matrixmod.h
matrixmod.c
```

This module owns:

```c
pmod_mat_print()
pmod_mat_fprint()

pmod_mat_zero()
pmod_mat_copy()
pmod_mat_identity()
pmod_mat_transpose()

pmod_mat_mul()
pmod_mat_mul_rect()

pmod_mat_vec_mul()
pmod_mat_transpose_vec_mul()

pmod_mat_set_col()
pmod_mat_get_col()

pmod_mat_linear_combination()
pmod_mat_diag_scale()
```

Do not add rank, kernel, span, or general pivot-search code to `matrixmod.c`.

## Square matrix multiplication

For new corank/CF code, prefer:

```c
pmod_mat_mul(C, A, B, n);
```

This computes:

```text
C = A B
```

for `n x n` square matrices.

The implementation uses a temporary output matrix and is expected to support exact output aliasing:

```c
pmod_mat_mul(A, A, B, n);
pmod_mat_mul(B, A, B, n);
pmod_mat_mul(A, A, A, n);
```

Do not assume support for arbitrary partially overlapping buffers.

## Rectangular multiplication

`pmod_mat_mul_rect()` remains for legacy code.

Its common inner dimension is:

```text
A_c == B_r
```

The multiplication loop must iterate over `A_c`, not `A_r`.

New Corank/CF code should not use the rectangular helper when all operands are square.

## Accumulation bounds

Several matrix helpers accumulate products in `uint64_t` and reduce modulo `MEDS_p` after a dot product or elementwise linear combination.

Before increasing dimensions or the modulus, verify bounds such as:

```text
n * (MEDS_p - 1)^2 <= UINT64_MAX
```

and, for a linear combination of `count` matrices:

```text
count * (MEDS_p - 1)^2 <= UINT64_MAX
```

Do not change parameters without reviewing all delayed-reduction bounds.

## Transpose

`pmod_mat_transpose()` supports:

- separate input and output;
- exact in-place transpose where `AT == A`.

Do not assume support for partial overlap.

## Matrix-vector multiplication

The matrix-vector functions use temporary output storage and are expected to support:

```c
pmod_mat_vec_mul(v, A, v, n);
pmod_mat_transpose_vec_mul(v, A, v, n);
```

Do not rely on other forms of partial overlap.

## Column helpers

Use:

```c
pmod_mat_set_col(A, col, v, n);
pmod_mat_get_col(v, A, col, n);
```

when building matrices whose columns are vectors.

Do not use `memcpy()` to insert a column into a row-major matrix.

Callers must guarantee:

```text
0 <= col < n
```

## Matrix linear combination

`pmod_mat_linear_combination()` computes:

```text
out = sum_i coeffs[i] * matrices[i]
```

where the matrices are stored consecutively.

For `n` slices:

```text
matrices[i * n * n + row * n + col]
```

This helper is intended for constructions such as:

```text
M(w) = sum_i w_i M_i
```

Callers must ensure that the total delayed accumulation fits in `uint64_t`.

Do not assume arbitrary partial-buffer overlap unless a dedicated test covers it.

## Diagonal scaling

`pmod_mat_diag_scale()` computes:

```text
out = diag(left_diag) * A * diag(right_diag)
```

elementwise in `O(n^2)` time.

Use it instead of constructing explicit diagonal matrices and performing two `O(n^3)` multiplications.

The implementation uses temporary output storage and is expected to support:

```c
pmod_mat_diag_scale(A, left, A, right, n);
```

---

# Pivot-aware elimination layer

## Ownership

Pivot search, reduced row-echelon form, rank, inverse, kernels, and span tests belong in:

```text
matrixelim.h
matrixelim.c
```

The module must contain only one shared pivot-aware elimination core.

Do not implement separate copies of Gaussian elimination for:

- rank;
- inverse;
- right kernel;
- left kernel;
- span tests;
- random invertible-matrix generation.

The internal RREF core should remain `static` unless a concrete external use is approved.

## Variable-time naming

The current elimination APIs use a `_vartime` suffix.

Examples:

```c
pmod_mat_rank_vartime()
pmod_mat_is_invertible_vartime()
pmod_mat_inv_vartime()
pmod_mat_right_kernel_corank1_vartime()
pmod_mat_left_kernel_corank1_vartime()
pmod_vec_in_span_vartime()
pmod_vec_extends_independent_set_vartime()
```

The suffix is semantically important.

These functions:

- search for pivots based on matrix entries;
- branch on zero and nonzero values;
- may skip zero factors;
- are not constant-time.

Do not remove `_vartime`.

Do not use these functions on secret-dependent inputs in a side-channel-sensitive path without a separate security review.

Current intended uses involve public or algorithmically revealed matrices in the Corank/CF computation. If that assumption changes, revisit the implementation.

## Rank

Use:

```c
int rank = pmod_mat_rank_vartime(A, n);
```

Expected range:

```text
0 <= rank <= n
```

The implementation must support skipped pivot columns.

Do not use legacy systematic-form elimination to compute rank.

## Invertibility

Use:

```c
pmod_mat_is_invertible_vartime(A, n)
```

Expected result:

```text
1 if rank(A) == n
0 otherwise
```

Do not compute a determinant solely to test invertibility.

## Inverse

Use:

```c
pmod_mat_inv_vartime(A_inv, A, n)
```

Expected result:

```text
0   success
-1  singular matrix
```

The implementation uses Gauss-Jordan elimination with an identity right-hand side.

Exact aliasing is expected to work:

```c
pmod_mat_inv_vartime(A, A, n);
```

because the input is copied before output is written.

Do not depend on the contents of the output buffer after a failed inversion unless the interface is explicitly strengthened.

## Corank-one right kernel

Use:

```c
pmod_mat_right_kernel_corank1_vartime(kernel, A, n)
```

It succeeds only when:

```text
rank(A) == n - 1
```

On success:

```text
A * kernel == 0
kernel != 0
```

It must reject:

- full-rank matrices;
- matrices of corank two or more;
- the zero matrix when `n > 1`.

The current implementation selects the first free column in the RREF pivot order and sets that free coordinate to one. This gives a deterministic representative.

Do not treat this as a general nullspace-basis function.

## Corank-one left kernel

Use:

```c
pmod_mat_left_kernel_corank1_vartime(kernel, A, n)
```

It succeeds only when:

```text
rank(A) == n - 1
```

On success:

```text
kernel^T * A == 0
kernel != 0
```

The current implementation obtains the left kernel by computing the right kernel of `A^T`.

Do not reverse left- and right-kernel calls in `BuildUVW`.

## Span membership

Use:

```c
pmod_vec_in_span_vartime(v, basis, basis_count, n)
```

to test:

```text
v in span(basis[0], ..., basis[basis_count - 1])
```

Callers must guarantee:

```text
0 <= basis_count <= n
```

The mathematical empty-span rule is:

```text
span(empty set) = { zero vector }
```

Do not pass invalid `basis_count` values and then depend on implementation fallback behavior.

## Extending an independent set

Use:

```c
pmod_vec_extends_independent_set_vartime(
    v, basis, basis_count, n);
```

only when the existing basis vectors are already known to be linearly independent.

Expected result:

```text
1  adding v increases rank from basis_count to basis_count + 1
0  v is dependent on the current set
```

Callers must guarantee:

```text
0 <= basis_count <= n
```

This helper is preferable inside `BuildUVW` when the algorithm maintains an independent list incrementally.

Do not use it as a replacement for the general `InSpan` predicate when the existing vectors may already be dependent.

---

# Legacy systematic-form helpers

The following functions are legacy MEDS helpers:

```c
pmod_mat_syst_ct()
pmod_mat_row_echelon_ct()
pmod_mat_back_substitution_ct()
```

They are intended for a matrix with a full-rank leading square block, such as:

```text
[A | B] -> [I | A^{-1} B]
```

They are not general rank or kernel routines.

Do not use them for:

- rank computation;
- corank-one detection;
- left kernel;
- right kernel;
- span membership;
- `BuildUVW`;
- `DiagonalNormalize`;
- `CF`.

The name `_ct` is historical. Do not assume these functions are formally constant-time solely because of the suffix.

Do not rename or delete them until all legacy call sites are audited.

When maintaining these functions:

- keep behavior unchanged unless a dedicated task says otherwise;
- eliminate `-Wshadow` warnings;
- use descriptive variable names such as `pivot_entry`, `target_row_entry`, and `reduced_entry`;
- avoid nested reuse of names such as `val`, `tmp0`, and `tmp1`;
- use an integer type wide enough for intermediate signed subtraction.

The compatibility wrapper:

```c
pmod_mat_inv(B, A, A_r, A_c)
```

may remain for legacy callers, but new square-matrix code should use:

```c
pmod_mat_inv_vartime(B, A, n)
```

---

# Random invertible matrices

`rnd_inv_matrix()` must reuse the shared matrix-elimination API.

Do not add another private Gaussian-elimination implementation to `util.c`.

Expected logic:

```text
sample random square matrix
while matrix is not invertible:
    sample again
```

using:

```c
pmod_mat_is_invertible_vartime()
```

For identical seeds and dimensions, generation must remain deterministic.

Any change to the acceptance test can change downstream KAT outputs. Treat such changes as algorithmically significant and document them.

---

# Trilinear-form preparation

The new algorithm represents a trilinear form through `n` square slices:

```text
M = (M_1, ..., M_n)
M_i in F_q^(n x n)
```

For a vector `w`:

```text
M(w) = sum_i w_i M_i
```

Use:

```c
pmod_mat_linear_combination()
```

to build `M(w)`.

For a vector `v`, the columns of `Phi_V(v)` are:

```text
M_i v
```

Use:

```c
pmod_mat_vec_mul()
pmod_mat_set_col()
```

For a vector `u`, the columns of `Phi_U(u)` are:

```text
M_i^T u
```

Use:

```c
pmod_mat_transpose_vec_mul()
pmod_mat_set_col()
```

Do not explicitly transpose every slice unless necessary.

When a new trilinear-form module is introduced, prefer dedicated files such as:

```text
triform.h
triform.c
```

Do not overload legacy `pi()` or `phi()` semantics without auditing all callers.

---

# Current triform-derived-operators phase

The current focused production phase is the standalone trilinear-form derived-operator layer.

This phase may add or modify only:

```text
triform.h
triform.c
Makefile        only to add triform.o
CMakeLists.txt  only to add triform.c
AGENTS.md       only to document the phase constraints
```

This phase must not modify:

```text
util.c legacy phi()
util.h legacy phi() declaration
meds.c phi() call sites
crypto_sign_keypair()
crypto_sign()
crypto_sign_open()
KeyGen, Sign, Verify mathematical flow
public-key, secret-key, or signature formats
params.h or api.h sizes
challenge parsing or hashing inputs
```

The trilinear form is stored as `n` consecutive row-major square slices:

```text
M[slice * n * n + row * n + col]
```

The public `triform` interfaces currently have the precondition:

```text
1 <= n && n <= MEDS_n
```

because the pivot-aware elimination layer still uses `MEDS_n`-bounded internal arrays.

`triform_matrix_at_w()` must compute:

```text
M(w) = sum_i w_i M_i
```

by calling `pmod_mat_linear_combination()`.

`triform_phi_u()` must build columns:

```text
column i = M_i^T u
```

using `pmod_mat_transpose_vec_mul()` followed by `pmod_mat_set_col()`.

`triform_phi_v()` must build columns:

```text
column i = M_i v
```

using `pmod_mat_vec_mul()` followed by `pmod_mat_set_col()`.

`triform_eval()` must satisfy:

```text
phi_M(u,v,w) = sum_i w_i u^T M_i v
             = u^T M(w) v
             = v^T Phi_U(u) w
             = u^T Phi_V(v) w
```

The corank-one wrappers must map directions exactly as follows and must not perform a separate rank computation before calling the kernel helper:

```text
Phi_U(u) -> pmod_mat_left_kernel_corank1_vartime()
Phi_V(v) -> pmod_mat_right_kernel_corank1_vartime()
M(w)     -> pmod_mat_left_kernel_corank1_vartime()
```

`triform_action_pullback()` is a new independent operator for:

```text
phi_out(x, y, z) = phi_M(Ax, By, Cz)
out_j = A^T * (sum_i C[i,j] M_i) * B
```

It must read `C` by columns for output slices:

```text
coefficient for input slice i and output slice j is C[i * n + j]
```

It must not call or wrap legacy `phi()`, and legacy `phi()` must not call it.

The public aliasing contract is:

```text
triform_matrix_at_w: out must not overlap M or w
triform_phi_u:       phi_u must not overlap M or u
triform_phi_v:       phi_v must not overlap M or v
triform_action_pullback:
                     out must not overlap M, A, B, or C
```

This phase does not implement `BuildUVW`, `DiagonalNormalize`, `CF`, or `Corank1Cal`, and it must not connect `triform_action_pullback()` to the current protocol.

---

# BuildUVW requirements

When `BuildUVW` is implemented, it should use the new pivot-aware helpers.

Expected primitive mapping:

```text
rank(Phi_U(u)) == n - 1
    -> pmod_mat_rank_vartime()
       or a combined corank-one kernel helper

LKer(Phi_U(u))
    -> pmod_mat_left_kernel_corank1_vartime()

rank(Phi_V(v)) == n - 1
    -> pmod_mat_rank_vartime()
       or a combined corank-one kernel helper

RKer(Phi_V(v))
    -> pmod_mat_right_kernel_corank1_vartime()

M(w)
    -> pmod_mat_linear_combination()

LKer(M(w))
    -> pmod_mat_left_kernel_corank1_vartime()

new vector independent of previous vectors
    -> pmod_vec_extends_independent_set_vartime()
```

Avoid computing the rank and then independently performing a second elimination for the kernel when a combined operation can be introduced cleanly.

If profiling shows repeated elimination is expensive, add a combined internal/public helper in `matrixelim`, not a duplicate implementation in the upper layer.

The vectors `u_i`, `v_i`, and `w_i` are columns of `U`, `V`, and `W`.

Use:

```c
pmod_mat_set_col()
```

to assemble them.

---

# DiagonalNormalize requirements

`DiagonalNormalize` should use:

```c
pmod_mat_mul()
pmod_mat_transpose()
pmod_mat_diag_scale()
GF_inv_checked()
GF_batch_inv()
```

Do not construct explicit diagonal matrices and multiply them as dense matrices.

Preferred normalization path:

1. Compute transformed slices.
2. Collect all required nonzero anchors.
3. Reject immediately if any anchor is zero.
4. Batch-invert independent anchor values where practical.
5. Apply row and column factors through `pmod_mat_diag_scale()`.
6. Apply per-slice scalar factors through field multiplication.

All anchor positions must come from the authoritative algorithm specification.

Do not guess or “simplify” anchor indices.

If the algorithm uses an index such as slice 5, validate that the parameter set satisfies the required minimum dimension.

---

# Canonical form requirements

`CF` is conceptually:

```text
CF = BuildUVW followed by DiagonalNormalize
```

The CF implementation must be deterministic for the same input form and starting point.

Determinism depends on:

- deterministic pivot selection;
- deterministic free-column selection;
- deterministic kernel representative;
- deterministic independence checks;
- deterministic normalization anchors;
- canonical serialization.

Do not introduce randomized pivot choices.

A successful CF implementation must be tested for:

- exact repeatability;
- invariance under the specified form action;
- invariance under allowed projective rescaling of the starting point;
- rejection when `BuildUVW` fails;
- rejection when a normalization anchor is zero.

---

# Legacy MEDS2endGen behavior

Until the upper-layer Corank/CF migration is explicitly performed, preserve the currently working KeyGen, Sign, Verify, serialization, and KAT behavior.

Do not silently alter:

- the existing `phi()` convention;
- multiplication order in responses;
- challenge parsing;
- key sizes;
- signature sizes;
- full-matrix hashing;
- public-key serialization;
- secret-key serialization.

If the latest Corank/CF specification conflicts with the current implementation, handle the migration in a dedicated plan and focused commits.

Do not mix a large upper-layer algorithm migration into a low-level linear-algebra fix.

---

# Serialization rules

Use existing bitstream conventions.

Every writer and reader must agree on:

- element count;
- element bit width;
- matrix boundary;
- `bs_finalize()` placement;
- row-major ordering.

Decoded field elements must satisfy:

```text
value < MEDS_p
```

Reject non-canonical encodings.

Do not confuse ordinary full field-element serialization with the removed MEDS public-key compression technique.

Whenever a serialized layout changes, update together:

```text
params.h
api.h
encoder
decoder
KAT expectations
```

---

# Constant-time and side-channel rules

Do not infer side-channel safety from a function name.

The following categories must be distinguished:

1. Explicitly reviewed constant-time operations.
2. Historical `_ct` functions not yet formally verified.
3. Explicit `_vartime` functions.

Any function with `_vartime` must be treated as variable-time.

Before using variable-time elimination on secret-dependent data:

- identify the secret;
- identify the observable branches and memory accesses;
- obtain explicit approval for the design;
- otherwise replace it with an appropriate constant-time method.

Do not remove security-significant suffixes during cosmetic refactoring.

---

# Error handling and contracts

Follow these conventions unless an existing interface explicitly differs:

```text
0   success
-1  failure
```

Boolean predicates use:

```text
1  true
0  false
```

Document preconditions in headers for:

- valid dimensions;
- nonzero inversion inputs;
- valid column indices;
- valid `basis_count`;
- non-overlapping buffers;
- allowed exact aliasing.

Do not silently accept invalid dimensions.

Do not use the content of an output buffer after a failed operation unless the interface explicitly guarantees its state.

---

# Coding style

Use C11-compatible code.

Prefer:

- descriptive local names;
- one logical operation per block;
- `const` on read-only pointers;
- explicit casts at type-width boundaries;
- `size_t` for byte counts;
- fixed row-major indexing;
- comments that explain mathematical intent.

Avoid:

- variable shadowing;
- generic names such as `val`, `tmp0`, and `tmp1` in nested scopes;
- duplicated modular arithmetic;
- duplicated Gaussian elimination;
- unrelated formatting changes;
- hard-coded parameter values;
- hidden changes to mathematical conventions.

The code should compile cleanly with:

```text
-Wall
-Wextra
-Wshadow
```

Where practical, also test:

```text
-Wpedantic
-Wconversion
-Wsign-conversion
```

Do not enable a warning as an error globally until existing unrelated warnings are understood, but new or modified code should not introduce warnings.

---

# Build rules

Inspect the build files before running commands.

The normal Make build is:

```sh
make clean
make
```

Run KAT:

```sh
make KAT
```

Run benchmark:

```sh
make BENCH
```

The current default build is expected to include:

```text
bench
PQCgenKAT_sign
```

Additional test executables may be conditionally present.

CMake build:

```sh
cmake -S . -B cmake-build-test -DCMAKE_BUILD_TYPE=Release
cmake --build cmake-build-test -j
```

Do not commit:

- `build/`;
- `debug/`;
- CMake build directories;
- object files;
- generated executables;
- temporary logs;
- generated KAT files unless explicitly requested.

---

# Testing strategy

## Production branch

The production branch should retain:

- production source code;
- benchmark code;
- KAT generator;
- essential build files.

Do not permanently add large functional, randomized, sanitizer, or stress-test suites to `master` unless explicitly requested.

## Test branch

Detailed tests should live on a separate branch, preferably:

```text
test/linear-algebra-primitives
```

or the current project test branch identified by the user.

Recommended test layout:

```text
tests/linear_algebra/
```

Test-only code should be separated from production source code.

## Fix branch

When a test reveals a production bug:

- keep the failing test on the test branch;
- fix production code on a dedicated `fix/...` branch;
- rebase the test branch onto the fix branch;
- do not mix the production fix and test harness into one commit.

## Required low-level tests

Finite-field tests:

- exhaustive add/sub/neg/mul where practical;
- inversion of every nonzero field element;
- checked rejection of zero inversion;
- batch inversion;
- canonical output range.

Matrix tests:

- zero/copy/identity;
- in-place and out-of-place transpose;
- square multiplication;
- rectangular multiplication regression;
- exact aliasing where documented;
- matrix-vector multiplication;
- transpose-matrix-vector multiplication;
- column insertion/extraction;
- linear combination;
- diagonal scaling.

Elimination tests:

- rank 0 through rank `n`;
- skipped pivot columns;
- invertibility;
- left and right inverse identities;
- strict corank-one kernel acceptance;
- rejection of full-rank and corank-two matrices;
- deterministic kernel representative;
- span membership;
- extension of an independent set.

Stability tests:

- fixed-seed randomized tests;
- ASan;
- UBSan;
- GCC and Clang;
- Debug and Release;
- repeated execution;
- KAT regression;
- API sign/verify smoke tests.

---

# Done criteria

A low-level field or matrix task is complete only when all relevant conditions hold:

- the project builds with Make;
- the project builds with CMake;
- no new compiler warnings are introduced;
- field outputs remain canonical;
- matrix storage remains row-major;
- square and rectangular multiplication pass reference tests;
- delayed-accumulation bounds remain valid;
- inverse satisfies both `A A^{-1} = I` and `A^{-1} A = I`;
- rank handles skipped pivot columns;
- corank-one kernels satisfy their defining equations;
- corank-one APIs reject other ranks;
- span and independence semantics match their documented contracts;
- documented exact aliasing works;
- ASan and UBSan pass;
- KAT behavior is unchanged unless an intentional algorithmic change is documented;
- production code and test-only code remain separated;
- no duplicate Gaussian elimination is introduced.

An upper-layer Corank/CF task is complete only when, in addition:

- `BuildUVW` is deterministic;
- all required corank-one checks are enforced;
- `U`, `V`, and `W` are assembled with vectors as columns;
- normalization rejects zero anchors;
- CF is deterministic;
- CF invariance tests pass;
- serialization is canonical;
- KeyGen, Sign, and Verify agree on exactly the same mathematical conventions and byte layout.

---

# Workflow rules for coding agents

For nontrivial work:

1. Inspect the current branch and latest commit.
2. Read this file.
3. Read the authoritative algorithm specification.
4. Inspect headers and all call sites before editing.
5. Produce a plan before a multi-file algorithm change.
6. Keep production changes and test-only changes separate.
7. Make focused commits.
8. Run targeted tests after each logical step.
9. Run full build and regression checks before completion.
10. Report assumptions and unresolved specification ambiguities.

Recommended commit separation:

```text
field arithmetic
ordinary matrix primitives
pivot-aware elimination
corank-one kernels
span helpers
legacy cleanup
tests
upper-layer BuildUVW
upper-layer DiagonalNormalize
CF integration
KeyGen/Sign/Verify migration
serialization and KAT updates
```

Do not:

- combine unrelated refactoring;
- change notation silently;
- change the meaning of `phi()` without auditing every call site;
- replace variable-time code with purported constant-time code without verification;
- add external dependencies without approval;
- commit generated artifacts;
- delete legacy helpers before all callers are migrated.

When reporting completion, summarize:

- branch used;
- base commit;
- files changed;
- interfaces added or modified;
- mathematical behavior changed;
- serialization or size changes;
- tests run;
- test results;
- remaining risks;
- specification assumptions.

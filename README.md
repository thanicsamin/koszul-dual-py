# Koszul Duality Utilities

This repository contains two command-line programs for formal calculations in
BGG category O:

- `category_o_projectives.py` computes projective covers and Verma flags.
- `koszul_duals.py` prints Koszul-dual block data and selected Koszul duality
  functor images.

The programs are intentionally explicit. They are useful for checking examples,
printing small tables, and keeping track of conventions. They do not construct
category O modules as Sage objects.

## Requirements

Both programs import SageMath through `sage.all`, so run them with a Python
environment where Sage is importable.

In this checkout, if Sage tries to write cache files under a read-only home
directory, point `DOT_SAGE` at a writable directory first:

```bash
mkdir -p .sage
export DOT_SAGE="$PWD/.sage"
```

Then run commands with `python3` or with the Sage-enabled Python in your local
virtual environment.

```bash
python3 koszul_duals.py --help
python3 category_o_projectives.py --help
```

The scripts support text output by default. Several modes also support
`--format json`, and `category_o_projectives.py` also supports `--format csv`.

## `koszul_duals.py`

Use this script when you want Koszul-dual data rather than just projective
covers. It has four main modes.

### Principal `sl2` Block

```bash
python3 koszul_duals.py sl2-principal
python3 koszul_duals.py sl2-principal --format json
```

This prints the explicit finite graded algebras from the notes:

- `A = End_O(P(-2) + P(0))`
- `B = Ext_A^*(L,L)`, where `L = L(0) + L(-2)`
- `E(B) = Ext_B^*(B_0,B_0)`

It also prints the vertex-swap self-duality isomorphism, projective resolutions,
and a small table of known functor images.

### `sl2` Functor Images

```bash
python3 koszul_duals.py sl2
python3 koszul_duals.py sl2 M0 M-2 L0 P-2 M0v
python3 koszul_duals.py sl2 'L(0)[0]' 'L(-2)[1]'
python3 koszul_duals.py sl2 'P(0)v[0]' 'I(0)[0]' 'I(0)v[0]'
python3 koszul_duals.py sl2 'M(0)v[0]' --format json
```

This computes images under the implemented functor
`K = RHom_O(-, L_block)`, where `L_block` is the sum of simples in the relevant
regular integral `sl2` block.

Accepted object labels include:

- `M(0)`, `M0`
- `M(-2)`, `M-2`
- `L(0)`, `L0`
- `P(-2)`, `P-2`
- `I(0)`, `I0`
- cohomological shifts such as `L(0)[1]`
- dominant dual Vermas such as `M(0)v`, `M(0)vee`, or `M(0)*`
- vee-duals of simples, projectives, and injectives: `L(0)v`, `P(0)v`,
  `I(0)v`

The vee suffix is normalized using the category O duality:

```text
L(lambda)v = L(lambda)
P(lambda)v = I(lambda)
I(lambda)v = P(lambda)
```

In regular `sl2`, the script also knows the injectives explicitly. The
dominant injective is the dominant dual Verma, and the antidominant injective is
projective-injective:

```text
I(n) = M(n)v
I(-n-2) = P(-n-2)
```

For shifts, the functor is contravariant, so the output shift changes sign.
For example:

```bash
python3 koszul_duals.py sl2 'L(0)[0]' 'L(-2)[1]'
```

prints:

```text
K(L(0)[0]) = P(-2)[0]
K(L(-2)[1]) = P(0)[-1]
```

The script does not parse `o+` or direct-sum notation as one argument. Compute
each summand separately and add the results. Thus:

```text
K(L(0)[0] o+ L(-2)[1]) = P(-2)[0] o+ P(0)[-1]
```

### Regular Finite-Type Blocks

```bash
python3 koszul_duals.py regular-block A2
python3 koszul_duals.py regular-block B3
python3 koszul_duals.py regular-block A2 --dual-map right
python3 koszul_duals.py regular-block A2 --format json
```

This reports the BGS Koszul-dual block for a regular integral finite type block.
The output is a projective-to-simple correspondence:

```text
P(w . lambda) -> L((w0 w) . lambda_dual)
```

by default. Pass `--dual-map right` to use `w w0` instead of `w0 w`.

For types `A`, `D`, `E`, `F4`, and `G2`, the Cartan type is self-dual. Types
`B` and `C` are interchanged.

This mode reports the dual block and vertex correspondence. It does not compute
the full higher-rank Ext multiplication table.

### `sl_n` Koszul Duality

```bash
python3 koszul_duals.py sln 3
python3 koszul_duals.py sln 3 1,0
python3 koszul_duals.py sl3
python3 koszul_duals.py sl3 0,1
python3 koszul_duals.py sl3 'P(0,1)[0]'
python3 koszul_duals.py sl3 'M(0,1)[0]'
python3 koszul_duals.py sl3 'L(-2,2)[1]'
python3 koszul_duals.py sl3 'I(0,1)v[0]'
```

For a bare weight, this prints the regular-block projective-to-simple
correspondence. Weights are written in fundamental-weight coordinates. For
example, in `sl3`, `1,0` means the weight `omega_1`.

For object inputs, this computes a concentrated functor image. Higher-rank
object mode currently supports `P(...)`, `L(...)`, `I(...)v`, and two special
Verma cases with optional shifts. Dominant Vermas are normalized as
`M(dominant) = P(dominant)`, and antidominant Vermas are normalized as
`M(antidominant) = L(antidominant)`. The input `I(...)v` is normalized to
`P(...)`.

The script does not currently support intermediate higher-rank Vermas or bare
injective images such as `I(...)` and `P(...)v`; those require
projective-resolution/Ext data that this script does not store.

As above, `--dual-map left` is the default and `--dual-map right` switches the
vertex correspondence.

## `category_o_projectives.py`

Use this script when you want formal projective covers, Verma flags, and
optional composition multiplicities.

### `sl2` Projectives

```bash
python3 category_o_projectives.py sl2 0
python3 category_o_projectives.py sl2 2 --composition-factors
python3 category_o_projectives.py sl2 -2
python3 category_o_projectives.py sl2 0 --all-in-block
python3 category_o_projectives.py sl2 0 --format json
```

For `sl2`, the dot action is:

```text
s . lambda = -lambda - 2
```

If `n >= 0`, the regular integral block has simples `L(n)` and `L(-n-2)`.

### `sl_n` Projectives

```bash
python3 category_o_projectives.py sln 3 0
python3 category_o_projectives.py sln 3 1,0 --composition-factors
python3 category_o_projectives.py sl3 1,0 --all-in-block
python3 category_o_projectives.py sl4 0,1,0
```

Weights are in fundamental-weight coordinates:

- `sl3 1,0` means `omega_1`
- `sl4 0,1,0` means `omega_2`
- `0` means the zero weight of the correct rank

This mode currently supports regular integral blocks. The output includes the
dot-action coordinate formula for each simple reflection.

### General Weyl Group Mode

```bash
python3 category_o_projectives.py weyl A2 --projective s1
python3 category_o_projectives.py weyl A3 --projective 1,2,1
python3 category_o_projectives.py weyl G2 --all --format json
python3 category_o_projectives.py weyl A2 --all --composition-factors
```

This mode uses Sage's Kazhdan-Lusztig polynomials to compute formal projective
covers in a regular integral block.

Projective labels are reduced words:

- `e` means the identity element.
- `1` means `s1`.
- `1,2,1` means `s1 s2 s1`.
- `w0` or `longest` means the longest Weyl group element.

The `--base-weight` option selects the indexing convention:

```bash
python3 category_o_projectives.py weyl A2 --projective s1 --base-weight dominant
python3 category_o_projectives.py weyl A2 --projective s1 --base-weight antidominant
```

The conventions are:

```text
dominant:     [Delta(x): L(w)] = P_{x,w}(1), x <= w
antidominant: [Delta(x): L(w)] = P_{w,x}(1), w <= x
```

## Common Workflows

Print the explicit `sl2` principal block algebra and its Koszul dual:

```bash
python3 koszul_duals.py sl2-principal
```

Find the image of a shifted simple:

```bash
python3 koszul_duals.py sl2 'L(0)[2]'
```

Find the image of a direct sum by listing its summands:

```bash
python3 koszul_duals.py sl2 'L(0)[0]' 'L(-2)[1]'
```

Compute all projective covers in an `sl2` block:

```bash
python3 category_o_projectives.py sl2 0 --all-in-block
```

Compute a regular `A2` block correspondence:

```bash
python3 koszul_duals.py regular-block A2
```

Compute formal projective covers for every element in a Weyl group:

```bash
python3 category_o_projectives.py weyl A2 --all
```

## Limitations

- SageMath must be importable before either script can run, because
  `koszul_duals.py` imports `category_o_projectives.py`.
- The `sl2` functor mode supports regular integral `sl2` blocks.
- Higher-rank `sl_n` object mode in `koszul_duals.py` supports `P(...)`,
  `L(...)`, `I(...)v` as an alias for `P(...)`, and `M(...)` only when the
  Verma is already projective or simple. It does not compute intermediate
  Verma images or bare injective images such as `I(...)` or `P(...)v`.
- Direct sums are handled by additivity. Pass the summands as separate object
  arguments.
- The higher-rank duality modes report block correspondences, not full Ext
  algebra multiplication tables.

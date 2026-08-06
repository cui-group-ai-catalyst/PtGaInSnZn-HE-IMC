# Generic formulas and bounded validation protocol

This note defines the reusable mathematical interface implemented by the P1
extension. It does not introduce a new thermodynamic equation. Its purpose is
to make the inputs, validation operations, outputs, and excluded inferences
explicit and machine-checkable.

## 1. Integer composition manifold

For `m` species sharing `s` integer sites, every composition vector satisfies

```text
n_1 + n_2 + ... + n_m = s,    n_i >= 0,
y_i = n_i / s.
```

The complete integer manifold contains

```text
N = C(s + m - 1, m - 1)
```

rows. The validator checks non-negative integer counts, a constant site sum,
unique vectors, required endmembers, and (when requested) complete-manifold
cardinality. Neither the element names nor the CSV row order enter this rule.

## 2. Pairwise CEF representation

For supplied per-atom energies `g_data(y)`, the implemented representation is

```text
g_CEF(y) = sum_i y_i g_i^0 + sum_(i<j) omega_ij y_i y_j,
```

where `g_i^0` is the supplied endmember energy and `omega_ij` is obtained by
least squares. For the manuscript Pt3X manifold, `omega_ij` is on the
per-alloy-atom energy basis. If a beta-sublattice parameter is used instead,
the manuscript convention is

```text
omega_ij^(atom) = Omega_ij^(beta-site) / 4.
```

This fit represents the supplied energy surface. It does not validate the
energy backend, establish the global ground state, or search competing phases.

## 3. Interpolation checks

For observations `g_k` and predictions `g_hat_k`, the reported error is

```text
RMSE = sqrt[(1/N) sum_k (g_hat_k - g_k)^2].
```

Three checks are generated from the same interface:

1. training residuals on all supplied rows;
2. an endmember-only ablation with all `omega_ij = 0`;
3. leave-one-composition-out validation for every non-endmember row.

If `group_holdout_element` is supplied, every integer count of that element is
held out as a complete composition family. This is a stricter internal
interpolation test, but it remains within the same host, prototype, elemental
space, structure-generation rule, and energy backend.

## 4. Energy-backend comparison

For every matched binary reference entry, three bundled columns are compared:
Materials Project DFT, UMA-s-1p1, and CHGNet. Each backend pair reports:

```text
RMSE, MAE, mean bias, Pearson r, and Spearman rho.
```

Spearman rho is the primary quantity because the screening use is candidate
ranking. Absolute RMSE is secondary: the reference set contains heterogeneous
M-Ga compounds and is not a direct validation set for the five-component
Pt3(Ga,In,Sn,Zn) landscape.

For a specified `k`, top-k agreement is reported as

```text
overlap count = |Top_k(A) intersect Top_k(B)|,
Jaccard = |intersection| / |union|.
```

For every pair of hosts `(p, q)`, a ranking reversal is counted when

```text
[E_A(p) - E_A(q)] [E_B(p) - E_B(q)] < 0.
```

The reversal fraction is the number of reversals divided by the number of
comparable, non-tied host pairs. The detailed CSV identifies each reversed
pair, so a reviewer can inspect the aggregate rather than accepting a single
correlation coefficient.

## 5. Claim boundary

The schema requires the following states for this release:

- software configurability: demonstrated for bundled inputs;
- interpolation: validated internally on the fixed L12 manifold;
- binary backend consistency: observed on the bundled reference set;
- new host or prototype transferability: not evaluated;
- N/O/S/P-containing compound transferability: not evaluated;
- synthesizability: not predicted.

Changing either of the final three states to a positive claim causes manifest
validation and the corresponding unit test to fail. A future scientific
extension must provide new structures, reference states, energy-method
validation, competing phases, interface/kinetic models, and independent
experimental evidence rather than only changing the JSON label.

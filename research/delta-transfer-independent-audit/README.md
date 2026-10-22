# Independent audit of the δ cross-cell halo transfer

**Conditional proof audit and exact finite controls; no new κ.** This answers
the request for a different contributor to check
[#4462](https://github.com/CrocSwap/integer-mult-bounds/pull/4462), at head
`89548dc186c9c90d459085b6d1be4834977705b9`, in
[handoffs #1166](https://github.com/CrocSwap/integer-mult-bounds/discussions/1166#discussioncomment-19111508).
The selected `23/10000` record and all predecessor packages are unchanged.

## Result and added coverage

The downward-transfer and bounded-excess conclusions hold under their stated
**row and column** hypotheses. The exponents **11, 24, 26, 26, 6** are sufficient
for the spectral endpoint, series length, bandwidth, agreement distance and
normalizer respectively. This audit does not claim their optimality.

The addition is an independently written proof/checker that:

- accounts for within-cell inflation **and** raw across-cell oracle error in
  one budget, rather than checking each separately against the full allowance;
- derives the nonnormal polynomial residual with the correct adjoint order;
- checks exact principal-section polynomial agreement and demonstrates that
  the inverse rows themselves are not exactly local;
- rejects row-only, missing-diagonal, wrong-adjoint and excess-budget controls;
- distinguishes the original weight regime from fine-prime setup precision.

No predecessor verifier is imported or copied. The generic polynomial identity,
sharp support radius and factor-two inverse comparison are credited to Nick
Petalas's [#4217 audit](../halo-locality-witness-audit/README.md), not new claims.

## Quantified contract

Let `p >= 2` be an integer, `0 <= W <= p-1`, and let `H` be a finite square
matrix with diagonal exactly one and physical bandwidth at most `W`. For cyclic
geometry, distance and restrictions must retain the physical periodic edges;
the executable fixtures below use noncyclic sections. Let `N` be entrywise
nonnegative off the diagonal. Assume **each** row and column off-sum of `N` is
at most `1-g`, where `g=1/(8p²)`.

Partition the retained off-diagonal entries into within-cell and across-cell
entries; all non-retained off-diagonal entries of `H` are zero. Set
`D=2^(32760p)`, `δ=4/D`, and `η=2^(-32768p)`. The downward route requires
`|H_ij| <= N_ij/(1-δ)` within cells and `|H_ij| <= N_ij` across cells.
The robust route instead requires the **total** row and column off-sums of
`|H|` to be at most `1-g+ρ`, with `0 <= ρ <= 1/(32p³)`.

These are hypotheses, not properties established here for all Gaussian callers.
Principal restrictions inherit unit diagonal, bandwidth and both off-sum bounds,
because restriction only deletes nonnegative terms from absolute sums.

## Proof

### 1. Downward mixed-cell transfer

For either a row or a column, write `T=T_in+T_out <= 1-g`. Then

```text
off(H) <= T_in/(1-δ) + T_out <= (1-g)/(1-δ) = s,
1-s = (g-δ)/(1-δ),
U = 1/(1-s) = 8p²(D-4)/(D-32p²).
```

The split uses no commutation or symmetry assumption. Both induced norms of
`H-I` are at most `s<1`, so the Neumann series gives both inverse norms at most
`U`, hence `||H^-1||₂ <= U`. Also `||H||₂ <= 1+s < 2`.

A short uniform ceiling is `8p² < U < 9p²`. To see this without constructing
huge rationals, `p² <= 4^p` for all integers `p>=2`, by induction. Therefore
`D >= 512·4^p >= 512p² > 288p²-32`. Clearing positive denominators proves
`U<9p²`; `32p²>4` proves the lower bound. Since `81 < 2^7 <= p^7` and
`9 < 2^4 <= p^4`, `U² < p^11` and `U < p^6`.

### 2. Robust excess and the combined raw-oracle bill

Write `q=g-ρ`. At the maximum allowed excess,

```text
q = (4p-1)/(32p³) >= 7/(64p²) > 1/(10p²).
```

Thus both inverse norms, and therefore the spectral inverse norm, are at most
`B=1/q < 10p²`. In particular `B² < 100p^4 < p^11`, `B < p^6`, and
`||H||₂ < 2`. The comparisons use `100<2^7` and `10<2^4`.

For a direct approximation to **all retained ideal entries**, a per-entry
absolute error at most `1/(64p^4)` suffices: each row or column has at most
`2W <= 2(p-1)` such entries, and their total is strictly below `1/(32p³)`.
This statement must not be read as an across-cell-only budget with the entire
allowance also spent on within-cell errors.

For the actual mixed construction, within-cell inflation contributes at most
`δ/(1-δ) <= 2δ`; a raw across-cell oracle with `|z-n|<2η` contributes at most
`4(p-1)η`. Their **sum** fits:

```text
2δ + 4(p-1)η = [8 + 4(p-1)/2^(8p)] / D < 9/D < 1/(32p³).
```

Here `4(p-1)<2^(8p)` follows by induction (successive left-side ratios are at
most two). Also `p³ <= 4^p`: the base is `8<=16`, and
`((p+1)/p)^3 <= 27/8 < 4`. Hence `D >= 512·4^p >= 512p³ > 288p³`.
This proves the combined budget for **all** integers `p>=2`, not just the finite
checked range. Downward repair still supplies the sharper ceiling `U`; the raw
oracle route establishes the robust ceiling instead.

### 3. Polynomial inverse and section agreement

Either route gives `||H||₂<2` and `||H^-1||₂<p^(11/2)`. Put

```text
α=p^-11, Q=I-α H*H, J=p^24,
P_J(H)=α sum_{j=0}^{J-1} Q^j H*.
```

The eigenvalues of `Q` lie in `[0,1-p^-22)`. Moreover
`α H*=(I-Q)H^-1`, so telescoping gives the exact, order-sensitive identity

```text
P_J(H)-H^-1 = -Q^J H^-1.
```

Both routes give `||H^-1||₂ < 10p²`, so
`||P_J(H)-H^-1||₂ < 10p² exp(-J p^-22) < p^6 exp(-p²)`.
The same statement holds for every principal section by the inherited bounds.
This is an analytic proof; the checker does not materialize `p^24` iterations.

Each term `Q^j H*` has bandwidth at most `(2j+1)W`. Thus the polynomial
bandwidth is at most `(2J-1)W`, and both that quantity and the conservative
path-agreement radius `(2J+1)W` are below `p^26`. For the latter,

```text
p^26-(2p^24+1)(p-1) = p^24[(p-1)^2+1]-(p-1) > 0.
```

A row farther than that path radius from the removed vertices has identical
global and embedded section **polynomial** rows. It follows by two remainder
bounds that the inverse-row difference is less than `2p^6 exp(-p²)`.
It is not zero in general. This distinction preserves the factor two already
identified in #4217; it does not change an exponent.

## Caller applicability and limitations

The source pins in [SOURCES.json](SOURCES.json) cover the original phase-cell
reference, LOG-A and FINE-PRIME. They authenticate bytes, not their mathematics.

- Original rad20 has `a_j>2^-8p` under its `u*d<=p` regime. Absolute error
  `<4η` then implies the stated `δ`. FINE-PRIME explicitly drops that weight
  regime: weights can be `exp(-O(p² log p))` and require
  `P_setup=Cp+64p² ceil(log₂p)`. The original `δ` must not silently be reused.
- FINE-PRIME's stated `poly(p) 2^{-(C-o(1))p}` row **and column** perturbation
  fits the robust budget eventually, conditional on the ideal gaps. This audit
  gives no explicit cutoff and does not derive those perturbation/gap premises.
- No page generator, Gaussian coefficient production, Schur propagation,
  tensor cost ledger, native controller or complete selected replay is checked.
- The five exponents are sufficient for this conditional route, not minimum
  exponents or a new multiplication theorem. No formal proof-kernel or
  independent human verification is claimed.

## Reproduce and evidence

The checker uses Python 3.12 and the standard library. It checks exact parameter
inequalities for `p=2..64`, the combined error budget, nonnormal matrix residuals,
section locality and meaningful failure controls. Its small matrix fixtures use
short series, not the astronomical series length in the quantified proof.

```sh
direnv allow
direnv exec . uv run --no-project python -B scripts/research_package.py \
  --package research/delta-transfer-independent-audit --check-only
direnv exec . uv run --no-project python -B scripts/research_package.py \
  --package research/delta-transfer-independent-audit
```

The root development environment pins Python 3.12; the same package can be run
with an existing Python 3.12 installation and no Nix dependency. Optimized
execution is also supported: all mathematical checks use explicit exceptions.

Local evidence on Python 3.12.14: metadata/pins passed; the registered replay
passed in ordinary and `-I -O -B` modes, each covering 63 parameter values,
5 positive fixtures and 11 adverse controls. [RESULT.json](RESULT.json) records
the observed direct optimized output; acceptance is recomputed, not read from
that receipt. These are agent-run results, not hosted CI or independent human
verification.

Orientation limitation: on restored canonical source bytes,
`scripts/research_baseline.py` passed source checking but then failed with
`KeyError: 'bit_lower'` when interpreting the newly selected productive-helper
certificate. This is separate from the focused audit; no selected cold replay
is claimed. An earlier inventory failure was caused by automatic formatting
during environment initialization; all of those unrelated edits were restored,
and the retained environment has no formatter or Git-hook configuration.

## Coordination and provenance

Started from canonical main `2a13eb1827784206202db90e5dc8327cf8a250db`, priority
revision 7. The explicit handoff overrides the generic #4100 task-card question;
[TASK.json](TASK.json) records this changed scope. Initial overlap refresh at
2026-10-22 00:38 UTC covered all 334 open PRs, 3,360 Discussion records and all
5,731 tracked top-level notes, with no coverage warnings. It found #4462,
#4473's uncovered-transfer receipt and the predecessor checks, not another
independent #4462 check. This is text coverage, not proof of novelty.

A scoped claim reply to the handoff was attempted through `post_update.py` and
the permitted GitHub transport; the mutation returned HTTP 403, and a read-back
showed no reply. No public reservation is claimed. Canonical main was fetched
again before submission and remained at the same commit.

The final overlap refresh at 2026-10-22 00:49 UTC covered 336 open PRs,
3,362 Discussion records and all 5,731 tracked top-level notes. It also found
James Reed's new [#4481](https://github.com/CrocSwap/integer-mult-bounds/pull/4481),
head `99df7cf237f7fdc801501afea3347f5ef9f74704`, which claims the next Schur-gap
step, not an independent check of #4462. That proposal is not audited here.

Pending sources are review targets, not runtime imports: #4462 at
`89548dc186c9c90d459085b6d1be4834977705b9`; #4452 at
`af68f5973a572305af174b3c49c440592172134d`; #4451 at
`be243a2818f0de3421c118dd21601aafe8e839dd`; #4449 at
`1fb77b377f570ca84a89a465b2ee4ce8674b7159`; #4442 at
`d079c2b058ebd3d7a9ff424113b095bf61175427`. No unmerged Git object is required
to replay this package.

## Open next questions

1. **Which actual page generator meets the joint error budget?** It matters
   because this conditional implication cannot certify stored coefficients.
   First pin a generator and measure or certify both row and column excess,
   including within-cell ratios, across-cell coefficients and diagonal handling.
   Blocker: no generator or all-size error proof is supplied by this audit.
2. **Does the proposed Schur continuation meet its physical caller?**
   [#4481](https://github.com/CrocSwap/integer-mult-bounds/pull/4481) now proposes
   the elimination-gap step. First independently audit that exact source,
   including Schur bandwidth, then its normalization/rounding and paid-layout
   bill. The core and complete cost remain unchanged here; this audit does not
   establish #4481's caller premises.

## Credit

Ryanson Jonathan (pinterbanget), #4462, supplies the audited transfer and
robustness claim and requests the independent check. James Reed (jwr404),
#4442/#4451/#4452, supplies the split, relative-error budget and exponent route.
Mihir S (mihir-s-05), #4449, supplies the across-cell rounding contract.
Nick Petalas, #4217, supplies the generic locality audit and factor-two
qualification. The retained rad20 authors/reviewers, Vlado Vojdanovski (#1673),
Chafik Boukhalfa, Cinna Davis and their source notices retain credit for the
Gaussian/fine-prime construction. David Harvey and Joris van der Hoeven retain
the underlying multiplication-source credit.

Prepared for **codyw912**, with substantial **OpenAI Codex** assistance, including
source scouting, a separately delegated checker and mathematical review.
New work is Apache-2.0. No predecessor file is changed.

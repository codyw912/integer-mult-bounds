#!/usr/bin/env python3
"""Independent finite exact checks; the all-p proof lives in README.md.

No predecessor implementation is imported, and no runtime files are written.
All mathematical acceptance tests use integers or Fractions, including under -O.
"""

import hashlib
import json
from fractions import Fraction as F
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def check_inputs():
    package = Path(__file__).resolve().parent
    root = package.parent.parent
    manifest = json.loads((package / "SOURCES.json").read_text())
    inputs = []
    for source in manifest["files"]:
        # Historical provenance belongs to the official package runner.
        if source["usage"] != "input":
            continue
        path = root / source["path"]
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        require(digest == source["sha256"], f"input SHA256 mismatch: {source['path']}")
        inputs.append(source["path"])
    require(bool(inputs), "source manifest has no live inputs")
    return inputs


def parameter_checks():
    for p in range(2, 65):
        D = 1 << (32760 * p)
        two8p = 1 << (8 * p)
        g, rho = F(1, 8 * p**2), F(1, 32 * p**3)
        # s = snum/common; 1-s = gapnum/common; U = common/gapnum.
        # Never normalize or square these million-bit rational denominators.
        common = 8 * p**2 * (D - 4)
        snum = (8 * p**2 - 1) * D
        gapnum = D - 32 * p**2
        require(D > 32 * p**2, f"p={p}: downward positive gap")
        require(snum + gapnum == common, f"p={p}: complementary split")
        require(8 * p**2 * gapnum < common, f"p={p}: U lower bound")
        require(D > 32 * (9 * p**2 - 1), f"p={p}: downward upper threshold")
        require(common < 9 * p**2 * gapnum, f"p={p}: U upper bound")
        require(81 < p**7 and 9 < p**4, f"p={p}: downward exponents 11,6")
        # delta/(1-delta) <= 2*delta when delta <= 1/2.
        require(D > 8, f"p={p}: within-cell inflation bound")
        # One shared row/column budget, not two independent rho budgets.
        require((8 * two8p + 4 * (p - 1)) * 32 * p**3 < D * two8p,
                f"p={p}: combined within/across raw-oracle excess")
        require(2 * (p - 1) * F(1, 64 * p**4) < rho,
                f"p={p}: per-entry approximation total")
        q = g - rho
        require(q == F(4 * p - 1, 32 * p**3), f"p={p}: robust gap")
        require(q >= F(7, 64 * p**2) > F(1, 10 * p**2),
                f"p={p}: robust gap lower bounds")
        B = 1 / q
        require(B < 10 * p**2 and B**2 < p**11 and B < p**6,
                f"p={p}: inverse norm and exponent bounds")
        require(0 < q < 1 and 1 + (1 - q) < 2,
                f"p={p}: both induced H norms below two")
        scale, J = F(1, p**11), p**24
        # Thus spec(Q) lies in [0, 1-p^-22); no numerical eigenvalues needed.
        require(4 * scale < 1 and scale * q**2 > F(1, p**22),
                f"p={p}: spectral contraction endpoints")
        require(0 < 1 - F(1, p**22) < 1, f"p={p}: exponential estimate domain")
        require(J * F(1, p**22) == p**2, f"p={p}: exponential error exponent")
        W = p - 1
        require((2 * J - 1) * W < p**26, f"p={p}: polynomial bandwidth")
        require((2 * J + 1) * W < p**26, f"p={p}: conservative agreement radius")
    return {
        "p_min": 2, "p_max": 64, "cases": 63, "arithmetic": "integer/Fraction",
        "downward": "8p^2 < U < 9p^2; U^2 < p^11; U < p^6",
        "combined_budget": "8/D + 4(p-1)/(D*2^(8p)) < 1/(32p^3)",
        "error": "B < p^6 and J*p^-22 = p^2; no astronomical polynomial evaluated",
        "bandwidth_and_agreement_radius": "(2J-1)W < p^26 and (2J+1)W < p^26",
    }


def identity(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def transpose(a):
    return [list(row) for row in zip(*a)]


def multiply(a, b):
    columns = transpose(b)
    return [[sum((x * y for x, y in zip(row, col)), F(0))
             for col in columns] for row in a]


def subtract(a, b):
    return [[x - y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scaled(a, scalar):
    return [[scalar * x for x in row] for row in a]


def power(a, exponent):
    result = identity(len(a))
    while exponent:
        if exponent % 2:
            result = multiply(result, a)
        exponent //= 2
        if exponent:
            a = multiply(a, a)
    return result


def inverse(a):
    n = len(a)
    rows = [row[:] + unit for row, unit in zip(a, identity(n))]
    for col in range(n):
        pivot = next((i for i in range(col, n) if rows[i][col]), None)
        if pivot is None:
            raise ValueError("singular matrix")
        rows[col], rows[pivot] = rows[pivot], rows[col]
        divisor = rows[col][col]
        rows[col] = [x / divisor for x in rows[col]]
        for i in range(n):
            if i != col:
                factor = rows[i][col]
                rows[i] = [x - factor * y for x, y in zip(rows[i], rows[col])]
    return [row[n:] for row in rows]


def off_norm(a):
    return max(sum(abs(x) for j, x in enumerate(row) if i != j)
               for i, row in enumerate(a))


def row_norm(a):
    return max(sum(abs(x) for x in row) for row in a)


def bandwidth(a):
    return max((abs(i - j) for i, row in enumerate(a)
                for j, x in enumerate(row) if x), default=0)


def admissible(a, p):
    return (all(a[i][i] == 1 for i in range(len(a)))
            and bandwidth(a) <= p - 1
            and max(off_norm(a), off_norm(transpose(a))) <= 1 - F(1, 8 * p**2))


def polynomial(a, J, scale):
    adjoint = transpose(a)  # All smoke fixtures are real rational matrices.
    Q = subtract(identity(len(a)), scaled(multiply(adjoint, a), scale))
    term, total = identity(len(a)), scaled(identity(len(a)), F(0))
    for _ in range(J):
        total = subtract(total, scaled(term, -1))
        term = multiply(term, Q)
    return scaled(multiply(total, adjoint), scale), Q, total


def matrix_checks():
    fixtures, controls = [], []
    a = [[F(1), F(1, 8), F(0)],
         [F(1, 16), F(1), F(1, 10)],
         [F(0), F(1, 20), F(1)]]
    at = transpose(a)
    require(admissible(a, 2), "3x3 unit diagonal, band and BOTH off norms")
    require(multiply(a, at) != multiply(at, a), "fixture must be nonnormal")
    inv = inverse(a)
    require(multiply(a, inv) == multiply(inv, a) == identity(3), "3x3 inverse")
    require(row_norm(inv) <= 1 / (1 - off_norm(a)), "inverse row bound")
    require(row_norm(transpose(inv)) <= 1 / (1 - off_norm(at)), "inverse column bound")
    fixtures.append("nonnormal_3x3_row_and_column_inverse_bounds")
    for J in (2, 3):
        P, Q, total = polynomial(a, J, F(1, 4))
        residual = scaled(multiply(power(Q, J), inv), -1)
        require(subtract(P, inv) == residual, f"J={J}: LEFT Q^J inverse identity")
        fixtures.append(f"nonnormal_residual_identity_J{J}")
        no_adjoint = scaled(multiply(total, a), F(1, 4))
        require(subtract(no_adjoint, inv) != residual, f"J={J}: missing adjoint control")
        wrong_Q = subtract(identity(3), scaled(multiply(a, at), F(1, 4)))
        wrong_total = subtract(identity(3), scaled(identity(3), 1))
        for j in range(J):
            wrong_total = subtract(wrong_total, scaled(power(wrong_Q, j), -1))
        wrong_P = scaled(multiply(wrong_total, at), F(1, 4))
        require(subtract(wrong_P, inv) != residual, f"J={J}: wrong Q orientation control")
        controls.extend([f"omitted_adjoint_J{J}", f"wrong_Q_orientation_J{J}"])

    n = 13
    h = identity(n)
    for i in range(n - 1):
        h[i][i + 1], h[i + 1][i] = F(1, 8), F(1, 16)
    require(admissible(h, 2), "locality fixture hypotheses")
    section = [row[1:12] for row in h[1:12]]
    P, _, _ = polynomial(h, 3, F(1, 4))
    Ps, _, _ = polynomial(section, 3, F(1, 4))
    embedded = [[F(0) for _ in range(n)] for _ in range(n)]
    for i in range(11):
        embedded[i + 1][1:12] = Ps[i]
    require(bandwidth(P) <= 5 and bandwidth(Ps) <= 5, "degree-five polynomial bandwidth")
    require(P[6] == embedded[6], "center polynomial row agrees including exterior zeros")
    require(P[1][1:12] != embedded[1][1:12], "near-edge INTERNAL polynomial row differs")
    invh, invs = inverse(h), inverse(section)
    require(multiply(h, invh) == multiply(invh, h) == identity(n), "full locality inverse")
    require(multiply(section, invs) == multiply(invs, section) == identity(11),
            "section locality inverse")
    require(invh[6][6] != invs[5][5], "exact inverse differs even at center diagonal")
    fixtures.append("principal_section_center_polynomial_locality")
    controls.extend(["near_edge_polynomial_disagreement", "inverse_not_finitely_local"])

    star = identity(18)
    for i in range(1, 18):
        star[i][0] = F(1, 2)
    p = 18  # Bandwidth is legal: the missing hypothesis is the column bound.
    require(bandwidth(star) <= p - 1 and off_norm(star) == F(1, 2), "row-only fixture")
    require(off_norm(transpose(star)) == F(17, 2) and not admissible(star, p),
            "row-only fixture must fail column hypothesis")
    # H=I+N, N^2=0; hence inverse=I-N, without a spectral solver.
    starinv = subtract(scaled(identity(18), 2), star)
    require(multiply(star, starinv) == multiply(starinv, star) == identity(18),
            "row-only exact inverse formula")
    require(sum(row[0]**2 for row in starinv) == F(21, 4) > 4,
            "row-only inverse spectral norm exceeds two")
    controls.append("row_only_hypothesis_insufficient")
    zero = [[F(0)]]
    require(off_norm(zero) == 0 and not admissible(zero, 2), "unit diagonal is required")
    try:
        inverse(zero)
    except ValueError as error:
        require(str(error) == "singular matrix", "unexpected inversion failure")
    else:
        raise ValueError("zero matrix incorrectly accepted as invertible")
    controls.append("missing_unit_diagonal")
    return fixtures, controls


def budget_controls():
    p, g, rho = 2, F(1, 32), F(1, 256)
    require((g - rho)**2 * p**11 > 1, "positive p=2 exponent-11 ceiling")
    require((g - F(1, 64))**2 * p**11 < 1, "overbudget exponent-11 ceiling failure")
    Tin = 1 - g
    delta = rho / (Tin + rho)
    inflation, across = Tin / (1 - delta) - Tin, rho
    require(inflation == rho and across == rho, "separate budgets each pass")
    require(inflation + across > rho, "combined separate-budget control must fail")
    D = 32 * p**2
    require(D - 32 * p**2 == 0, "downward gap boundary is singular")
    D = 64 * p**2
    U = F((D - 4) * 8 * p**2, D - 32 * p**2)
    require(U > 9 * p**2, "insufficient D must fail downward upper estimate")
    return ["p2_positive_robust_budget"], [
        "overbudget_exponent11_ceiling_not_gaussian_counterexample",
        "separate_budgets_pass_but_combined_budget_fails",
        "downward_zero_gap_boundary", "downward_insufficient_D_upper_bound",
    ]


def main():
    inputs = check_inputs()
    parameters = parameter_checks()
    fixtures, controls = matrix_checks()
    budget_fixtures, budget_failures = budget_controls()
    fixtures.extend(budget_fixtures)
    controls.extend(budget_failures)
    print(json.dumps({
        "status": "PASS_DELTA_TRANSFER_INDEPENDENT_AUDIT",
        "claim": "conditional PR4462 at 89548dc186c9c90d459085b6d1be4834977705b9",
        "parameters": parameters,
        "fixtures": {"count": len(fixtures), "names": fixtures},
        "negative_controls": {"count": len(controls), "names": controls},
        "combined_budget": "checked jointly for every integer p=2..64",
        "input_sha256": {"count": len(inputs), "paths": inputs},
        "no_new_kappa": True,
        "limitations": [
            "Finite p=2..64 checks are not the all-p proof in README.md.",
            "Conditional matrix hypotheses, not the full Gaussian supplier or multiplication theorem.",
            "Small rational real matrices test identities and locality, not astronomical J=p^24.",
            "Overbudget control rejects this exponent route; it is not a Gaussian counterexample.",
            "Combined-budget adverse fixture is arithmetic, not a realized Gaussian matrix.",
            "Historical provenance is checked by the official runner, not live byte comparison.",
        ],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

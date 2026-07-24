"""Reproduce the frozen ZHP7C Bayesian prior catalog for clean-t0 v2.

This is a normal summary-likelihood analysis. It propagates uncertainty from
the reported HR/RMST standard errors, but not the full estimation uncertainty
from fitting the CCW-IPCW models.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from statistics import NormalDist


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
REPORT = RESULTS / "strictK3_v2_bayes_report.md"
TSV = RESULTS / "strictK3_v2_bayes_prior_catalog_results.tsv"

NORMAL = NormalDist()


def cdf(value: float, mean: float, sd: float) -> float:
    return NormalDist(mean, sd).cdf(value)


def conjugate_update(estimate: float, se: float, prior_mean: float, prior_sd: float):
    variance = 1.0 / (1.0 / se**2 + 1.0 / prior_sd**2)
    mean = variance * (estimate / se**2 + prior_mean / prior_sd**2)
    return mean, math.sqrt(variance)


HR_PRIORS = [
    ("Data-dominant", 0.000, 1.000, "primary"),
    ("Neutral", 0.000, 0.355, "primary"),
    ("Skeptical", 0.000, 0.205, "primary"),
    ("Mildly optimistic", -0.105, 0.205, "sensitivity"),
    ("Mildly pessimistic", 0.095, 0.205, "sensitivity"),
    ("Discounted-historical neutral", 0.000, 0.250, "sensitivity"),
    ("ELAIN-informed stress test", -0.223, 0.300, "sensitivity"),
]

RMST_PRIORS = [
    ("Data-dominant", 0.000, 3.000, "primary"),
    ("Neutral/skeptical", 0.000, 1.000, "primary"),
    ("Mildly optimistic", 1.000, 2.000, "sensitivity"),
    ("Mildly pessimistic", -1.000, 2.000, "sensitivity"),
    ("Discounted-historical neutral", 0.000, 0.750, "sensitivity"),
]


def hr_results():
    estimate = math.log(1.044)
    se = (math.log(1.181) - math.log(0.922)) / (2 * 1.959964)
    rows = []
    for name, prior_mean, prior_sd, role in HR_PRIORS:
        mean, sd = conjugate_update(estimate, se, prior_mean, prior_sd)
        benefit = cdf(math.log(0.900), mean, sd)
        equivalence = cdf(math.log(1.100), mean, sd) - benefit
        harm = 1.0 - cdf(math.log(1.100), mean, sd)
        rows.append(
            {
                "endpoint": "HR",
                "prior": name,
                "role": role,
                "prior_mean": prior_mean,
                "prior_sd": prior_sd,
                "posterior_location": math.exp(mean),
                "posterior_sd": sd,
                "benefit": benefit,
                "equivalence": equivalence,
                "harm": harm,
            }
        )
    return estimate, se, rows


def rmst_results():
    estimate = -0.66
    se = (1.44 - (-2.76)) / (2 * 1.959964)
    rows = []
    for name, prior_mean, prior_sd, role in RMST_PRIORS:
        mean, sd = conjugate_update(estimate, se, prior_mean, prior_sd)
        benefit = 1.0 - cdf(1.0, mean, sd)
        equivalence = cdf(1.0, mean, sd) - cdf(-1.0, mean, sd)
        harm = cdf(-1.0, mean, sd)
        rows.append(
            {
                "endpoint": "RMST",
                "prior": name,
                "role": role,
                "prior_mean": prior_mean,
                "prior_sd": prior_sd,
                "posterior_location": mean,
                "posterior_sd": sd,
                "benefit": benefit,
                "equivalence": equivalence,
                "harm": harm,
            }
        )
    return estimate, se, rows


def percent(value: float) -> str:
    return f"{100 * value:.1f}%"


def write_outputs() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    hr_estimate, hr_se, hr_rows = hr_results()
    rmst_estimate, rmst_se, rmst_rows = rmst_results()
    rows = hr_rows + rmst_rows

    with TSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# Strict physiologic-KDIGO3 v2 Bayesian summary-likelihood reanalysis",
        "",
        "## Registration and execution boundary",
        "",
        "The MIMIC-IV effect estimate was completed before OSF BCD52. The Bayesian reanalysis was the sole prospectively registered component and was executed without further database access. The prior parameters below follow the frozen ZHP7C catalog. MAP and robust-mixture priors were not executed because their numerical construction was not frozen; this deviation is disclosed rather than back-filled post hoc.",
        "",
        "## Summary likelihoods",
        "",
        f"- log-HR estimate: {hr_estimate:.5f}; SE: {hr_se:.5f}",
        f"- RMST difference: {rmst_estimate:.2f} days; SE: {rmst_se:.4f}",
        "",
        "## HR posterior regions",
        "",
        "Important benefit was HR<0.900, practical equivalence HR 0.900-1.100, and important harm HR>1.100.",
        "",
        "| Prior | Role | Posterior median HR | Important benefit | Practical equivalence | Important harm |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in hr_rows:
        lines.append(
            f"| {row['prior']} | {row['role']} | {row['posterior_location']:.3f} | "
            f"{percent(row['benefit'])} | {percent(row['equivalence'])} | {percent(row['harm'])} |"
        )

    lines.extend(
        [
            "",
            "## RMST posterior regions",
            "",
            "Important benefit was an RMST difference >=1 day, practical equivalence -1 to 1 day, and important harm <=-1 day.",
            "",
            "| Prior | Role | Posterior mean, days | Important benefit | Practical equivalence | Important harm |",
            "|---|---|---:|---:|---:|---:|",
        ]
    )
    for row in rmst_rows:
        lines.append(
            f"| {row['prior']} | {row['role']} | {row['posterior_location']:.3f} | "
            f"{percent(row['benefit'])} | {percent(row['equivalence'])} | {percent(row['harm'])} |"
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "Across the three primary HR priors, the probability of clinically important benefit was 0.8%-0.9%, practical-equivalence probability was 78.8%-81.5%, and clinically important harm probability was 17.7%-20.3%. The mildly optimistic and ELAIN-informed stress-test priors yielded only 1.2% and 1.3% probability of important benefit, respectively. These posterior probabilities are conditional on the IPCW-adjusted summary likelihood and do not remove residual confounding or fully propagate weight-model uncertainty.",
            "",
        ]
    )
    REPORT.write_text("\n".join(lines), encoding="utf-8")

    for row in rows:
        location = row["posterior_location"]
        print(
            f"{row['endpoint']}\t{row['prior']}\t{location:.3f}\t"
            f"{percent(row['benefit'])}\t{percent(row['equivalence'])}\t{percent(row['harm'])}"
        )
    print(f"REPORT={REPORT}")
    print(f"TSV={TSV}")


if __name__ == "__main__":
    write_outputs()

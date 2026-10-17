import csv
import os
with open(os.environ["CH01_PRIMARY_CSV"]) as handle:
    input_rows = [r for r in csv.DictReader(handle) if r["population"] == "Equipoise SOFA 8-14"]
assert len(input_rows) == 1
CURRENT_HR = float(input_rows[0]["hr"])
CURRENT_CI_LOWER = float(input_rows[0]["ci_lower"])
CURRENT_CI_UPPER = float(input_rows[0]["ci_upper"])
import math
from statistics import NormalDist

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

def hr_results():
    estimate = math.log(CURRENT_HR)
    se = (math.log(CURRENT_CI_UPPER) - math.log(CURRENT_CI_LOWER)) / (2 * 1.959964)
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

if __name__ == "__main__":
    import json
    estimate, se, results = hr_results()
    print(json.dumps({"log_hr": estimate, "se": se, "rows": results}, indent=2))

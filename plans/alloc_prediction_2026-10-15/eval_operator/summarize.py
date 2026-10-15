"""Summarize saved predictions and decisions without refitting a model."""

import csv
from pathlib import Path
from statistics import mean, median

out = Path(__file__).resolve().parent


def read(name):
    with (out / name).open() as stream:
        return list(csv.DictReader(stream))


lines = ["# Computed results", "", "These are offline predictions, not measured speedups.", ""]
r = read("prediction.csv")
lines += [
    "## Prediction from 100 draws per chain",
    "",
    "Tau ratios use only the following 300 draws and the eight slowest coordinates in those draws. "
    "Lag-10 errors use the immediately following 100 draws. No evaluation draw enters the training covariance.",
    "",
    "| Dataset | Predictor | Median lag-10 absolute error | Median predicted / future tau, slowest eight |",
    "|---|---|---:|---:|",
]
for dataset in ("PM burn dirs", "pair_new_2 dirs", "pair_new_3 dirs", "pair_point_2 dirs", "pair_point_3 dirs", "full_T30 sample", "full_forst sample"):
    for method in ("scalar", "joint_group_0", "joint_aug_0", "joint_quad5"):
        x = [row for row in r if row["dataset"] == dataset and row["method"] == method]
        if not x:
            continue
        future = [row for row in x if row["future_draws"] == "300" and row["slow"] == "TRUE"]
        lag10 = median(abs(float(row["pred10"]) - float(row["real10"])) for row in x)
        ratio = median(float(row["pred_tau"]) / float(row["run_tau"]) for row in future)
        lines.append(f"| {dataset} | {method} | {lag10:.3f} | {ratio:.3f} |")
lines += [
    "",
    "## First three PM burn decisions",
    "",
    "These use the allocator's mean directions fixed from the checkpoint, with the same basis for every chain.",
    "",
    "| Training block | Predictor | Slowest-eight predicted / future-300 tau | Lag-10 absolute error |",
    "|---:|---|---:|---:|",
]
for b in ("1", "2", "3"):
    for method in ("scalar", "joint_group_0", "joint_aug_0", "joint_quad5"):
        x = [row for row in r if row["dataset"] == "PM burn dirs" and row["method"] == method and row["block"] == b]
        slow = [row for row in x if row["slow"] == "TRUE"]
        ratio = median(float(row["pred_tau"]) / float(row["run_tau"]) for row in slow)
        mae = median(abs(float(row["pred10"]) - float(row["real10"])) for row in x)
        lines.append(f"| {b} | {method} | {ratio:.3f} | {mae:.3f} |")
replays = []
for name in ("replay.csv", "replay_subject.csv", "replay_conditional.csv", "replay_latest.csv", "replay_checkpoint.csv"):
    if (out / name).exists():
        replays.extend(read(name))
lines += [
    "",
    "## Decisions from identical saved measurements",
    "",
    "Costs are predicted reserved CPU seconds per iteration, averaged over chains. "
    "Particle rates are the shared maximum over chains. HMC rates include regime rates.",
    "",
    "| State | Variant | Cost | Particle rate | HMC rate | Shear rounds |",
    "|---|---|---:|---:|---:|---:|",
]
for state in dict.fromkeys(row["state"] for row in replays):
    for variant in ("baseline", "point", "point_no_credit", "kappa", "kappa_no_credit", "own_subject", "baseline_subject", "conditional", "latest"):
        x = [row for row in replays if row["state"] == state and row["variant"] == variant]
        if x:
            vals = [mean(float(row[k]) for row in x) for k in ("cost", "q", "hmc", "shear")]
            lines.append(f"| {state} | {variant} | " + " | ".join(f"{v:.3f}" for v in vals) + " |")
(out / "RESULTS.md").write_text("\n".join(lines) + "\n")
print(out / "RESULTS.md")

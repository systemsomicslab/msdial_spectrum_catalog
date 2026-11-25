"""Draft figures for the WS8 P3 diagnostic-rule measurement (not the publication figures).

Inputs:
- ``--summary``: the unsuppressed per-rule summary of a frozen ``p3_measure.py`` run
  (``rules_summary.local-vs21.tsv``: pooled / vs20 / vs21 cluster counts, synthetic precision,
  per instrument class). The run's own CSVs cannot be used for this, because their held-out
  suppression also blanks development-set and passing-rule cells (a known issue).
- ``--msfinder-run``: optionally, a ``p3_msfinder.py`` run folder (development set).

Output: ``p3_summary.local-vs21.{png,svg}`` and ``p3_msfinder_offon.tsv``. Panels:

A. Cluster-level precision with Wilson 95 % intervals for every measured rule that fired.
B. The rules with a passing label: precision by source, on the synthetic sets and by instrument class.
C. MS-FINDER structure ranking, rules off against each weight.
D. Gate labels.

Everything that reads held-out values is named ``*.local-vs21.*`` and must never be committed.
Held-out cells with fewer than 5 fired clusters are suppressed in the figure.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

SUPPRESS_BELOW = 5
PASSING = ("pass_strict", "pass_provisional_only")
LABEL_COLORS = {
    "pass_strict": "#1b7837",
    "pass_provisional_only": "#5aae61",
    "incomplete_provisional": "#c2a5cf",
    "fail": "#c51b7d",
    "not_estimable": "#9e9e9e",
    "not_measured": "#d9d9d9",
}
INSTRUMENTS = ("Orbitrap/FT", "QTOF-type", "other")
TRUTHY = {"true", "1", "1.0"}


def wilson(k: float, n: float, z: float = 1.959963984540054) -> tuple[float, float]:
    """Wilson score interval for k successes out of n."""
    if not n or math.isnan(n):
        return float("nan"), float("nan")
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return centre - half, centre + half


def short_rule(rule_id: str) -> str:
    family, name, mode = rule_id.split(".")
    return f"{name if name == family else family + '/' + name} ({mode})"


def _flag(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin(TRUTHY)


def load_summary(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    for col in frame.columns:
        if col.endswith(("_c_fires", "_c_tp", "_c_fp", "_c_precision", "_c_wilson_lo")) or col == "prior_stress":
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
    return frame


def msfinder_summary(per_spectrum: Path) -> pd.DataFrame:
    t = pd.read_csv(per_spectrum, sep="\t", dtype=str)
    t = t[t["set"] == "selected"]
    rows = []
    for (family, variant), g in t.groupby(["family", "variant"]):
        defined = _flag(g["class_defined"])
        rows.append({
            "family": family,
            "variant": variant,
            "n": len(g),
            "top1": _flag(g["top1"]).mean(),
            "top5": _flag(g["top5"]).mean(),
            "class_acc": _flag(g.loc[defined, "class_correct"]).mean() if defined.any() else float("nan"),
            "boost_wrong": _flag(g["boost_wrong"]).mean(),
        })
    return pd.DataFrame(rows)


def _point(ax, y: float, tp: float, fires: float, color: str, marker: str = "o") -> None:
    lo, hi = wilson(tp, fires)
    p = tp / fires
    ax.errorbar([p], [y], xerr=[[max(0.0, p - lo)], [max(0.0, hi - p)]], fmt=marker, color=color, ecolor=color, capsize=2, ms=4)
    ax.text(1.01, y, f"{int(tp)}/{int(fires)}", va="center", fontsize=7)


def _thresholds(ax) -> None:
    for x, style in ((0.90, ":"), (0.95, "--")):
        ax.axvline(x, color="#555555", ls=style, lw=0.8)


def panel_fired(ax, summary: pd.DataFrame) -> None:
    fired = summary[summary["pooled_c_fires"] > 0].copy()
    fired["p"] = fired["pooled_c_tp"] / fired["pooled_c_fires"]
    fired = fired.sort_values(["p", "pooled_c_fires"]).reset_index(drop=True)
    for i, row in fired.iterrows():
        _point(ax, i, row["pooled_c_tp"], row["pooled_c_fires"], LABEL_COLORS.get(row["label"], "#000000"))
    ax.set_yticks(range(len(fired)))
    ax.set_yticklabels([short_rule(r) for r in fired["rule_id"]], fontsize=7)
    _thresholds(ax)
    ax.set_xlim(0, 1.08)
    ax.set_xlabel("cluster-level precision with Wilson 95 % (development + held-out pooled); colour = gate label")
    ax.set_title(f"A. Measured rules that fired ({len(fired)} of {len(summary)}); "
                 "dotted 0.90 provisional, dashed 0.95 strict", fontsize=9)


def panel_passing(ax, summary: pd.DataFrame) -> None:
    passing = summary[summary["label"].isin(PASSING)].sort_values("rule_id")
    labels, y = [], 0.0
    for _, row in passing.iterrows():
        name = short_rule(row["rule_id"])
        for src, text in (("pooled", "pooled"), ("vs20", "development"), ("vs21", "held-out")):
            fires, tp = row[f"{src}_c_fires"], row[f"{src}_c_tp"]
            if src == "vs21" and (pd.isna(fires) or fires < SUPPRESS_BELOW):
                ax.text(0.02, y, f"n < {SUPPRESS_BELOW}, suppressed", va="center", fontsize=7, color="#777777")
            elif fires > 0:
                _point(ax, y, tp, fires, "#2166ac")
            labels.append((y, f"{name} | {text}"))
            y += 1
        for col, text in (("chimeric_0.1_c_precision", "chimeric 10 %"), ("chimeric_0.3_c_precision", "chimeric 30 %"),
                          ("chimeric_0.5_c_precision", "chimeric 50 %"), ("insource_c_precision", "in-source recast")):
            value = pd.to_numeric(row[col], errors="coerce")
            if pd.notna(value):
                ax.plot([value], [y], "D", color="#b2182b", ms=4)
                ax.text(1.01, y, f"{value:.3f}", va="center", fontsize=7)
            labels.append((y, f"{name} | {text}"))
            y += 1
        for inst in INSTRUMENTS:
            p = pd.to_numeric(row[f"{inst}_c_precision"], errors="coerce")
            n = pd.to_numeric(row[f"{inst}_c_fires"], errors="coerce")
            if pd.notna(p) and pd.notna(n) and n > 0:
                _point(ax, y, round(p * n), n, "#762a83", "s")
            labels.append((y, f"{name} | {inst}"))
            y += 1
        y += 0.6
    ax.set_yticks([pos for pos, _ in labels])
    ax.set_yticklabels([text for _, text in labels], fontsize=7)
    ax.invert_yaxis()
    _thresholds(ax)
    ax.set_xlim(0, 1.1)
    ax.set_xlabel("cluster-level precision\nblue: source (Wilson 95 %); red: synthetic; purple: instrument class")
    ax.set_title("B. Rules with a passing label (author-reviewed in rule set 0.2.0)", fontsize=9)


def panel_msfinder(ax, ms: pd.DataFrame, run_name: str) -> None:
    variants = sorted(ms["variant"].unique(), key=lambda v: (v != "off", v))
    metrics = (("top1", "top-1"), ("top5", "top-5"), ("class_acc", "ChemOnt class"))
    families = sorted(ms["family"].unique())
    width = 0.8 / len(variants)
    ticks, names = [], []
    for fi, family in enumerate(families):
        n = int(ms.loc[ms["family"] == family, "n"].iloc[0])
        for mi, (col, metric) in enumerate(metrics):
            base = fi * (len(metrics) + 1) + mi
            for vi, variant in enumerate(variants):
                value = ms.loc[(ms["family"] == family) & (ms["variant"] == variant), col].iloc[0]
                shade = 0.3 + 0.65 * vi / max(1, len(variants) - 1)
                ax.bar(base + vi * width, value, width, color=plt.cm.Blues(shade),
                       label=variant if (fi, mi) == (0, 0) else None)
                ax.text(base + vi * width, value + 0.01, f"{value:.2f}", ha="center", fontsize=5, rotation=90)
            ticks.append(base + 0.4 - width / 2)
            names.append(metric)
        ax.text(fi * (len(metrics) + 1) + 1 + 0.4 - width / 2, -0.09, f"{family} (n={n})", ha="center",
                fontsize=8, transform=ax.get_xaxis_transform())
    ax.set_xticks(ticks)
    ax.set_xticklabels(names, fontsize=7)
    ax.set_ylim(0, 1.1)
    ax.set_ylabel("fraction of selected spectra")
    ax.legend(fontsize=7, title="ranking", title_fontsize=7, ncol=len(variants))
    ax.set_title(f"C. MS-FINDER structure ranking, development set ({run_name})", fontsize=9)


def panel_labels(ax, summary: pd.DataFrame) -> None:
    counts = summary["label"].value_counts()
    order = [k for k in LABEL_COLORS if k in counts.index]
    ax.barh(order, [counts[k] for k in order], color=[LABEL_COLORS[k] for k in order])
    for i, k in enumerate(order):
        ax.text(counts[k] + 0.3, i, str(counts[k]), va="center", fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("rules")
    ax.set_title(f"D. Gate labels ({len(summary)} measured rules)", fontsize=9)


def draw(summary_path: Path, msfinder_run: Path | None, out: Path, title: str) -> None:
    summary = load_summary(summary_path)
    ms = msfinder_summary(msfinder_run / "per_spectrum.tsv") if msfinder_run else pd.DataFrame()
    fig, axes = plt.subplots(2, 2, figsize=(16, 12), gridspec_kw={"width_ratios": [1.15, 1]})
    fig.suptitle(title, fontsize=12)
    panel_fired(axes[0, 0], summary)
    panel_passing(axes[0, 1], summary)
    if ms.empty:
        axes[1, 0].set_axis_off()
    else:
        panel_msfinder(axes[1, 0], ms, msfinder_run.name)
    panel_labels(axes[1, 1], summary)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out.mkdir(parents=True, exist_ok=True)
    stem = "p3_summary.local-vs21"
    fig.savefig(out / f"{stem}.png", dpi=150)
    fig.savefig(out / f"{stem}.svg")
    plt.close(fig)
    if not ms.empty:
        ms.to_csv(out / "p3_msfinder_offon.tsv", sep="\t", index=False)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--summary", type=Path, required=True, help="unsuppressed rules_summary.local-vs21.tsv")
    parser.add_argument("--msfinder-run", type=Path, help="p3_msfinder run folder (development set)")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--title", default="WS8 P3 diagnostic-ion rules: draft summary")
    args = parser.parse_args(argv)
    draw(args.summary, args.msfinder_run, args.out, args.title)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Check BPM example output; optionally plot it with matplotlib."""
import argparse
import csv
import json
import math
from pathlib import Path

FREQUENCIES = (100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000,
               1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000,
               10000, 12500, 16000, 20000, 25000, 31500, 40000)
ACTIVE = ("LBL_dB", "TBL_pressure_dB", "TBL_suction_dB", "TBL_separation_dB")
INACTIVE = ("bluntness_dB", "tip_dB", "inflow_dB")


def analyze(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(handle)]
    if tuple(row["frequency_Hz"] for row in rows) != FREQUENCIES:
        raise ValueError("Expected the 27 Table D2 band centers in order")
    for row in rows:
        for key in (*ACTIVE, "total_dB", "reference_D2_dB", "difference_dB"):
            if not math.isfinite(row[key]):
                raise ValueError(f"Non-finite {key} at {row['frequency_Hz']} Hz")
        if any(row[key] != -math.inf for key in INACTIVE):
            raise ValueError("Disabled mechanisms must have zero energy (-inf dB)")
        total = 10 * math.log10(sum(10 ** (row[key] / 10) for key in ACTIVE))
        if abs(total - row["total_dB"]) > 1e-8:
            raise ValueError("Total SPL does not equal the sum of source energies")
        if abs(row["total_dB"] - row["reference_D2_dB"] - row["difference_dB"]) > 1e-8:
            raise ValueError("Incorrect reported difference")
    errors = [row["difference_dB"] for row in rows]
    peak = max(rows, key=lambda row: row["total_dB"])
    worst = max(rows, key=lambda row: abs(row["difference_dB"]))
    if abs(worst["difference_dB"]) > .01:
        raise ValueError("Table D2 comparison exceeds 0.01 dB")
    return rows, {
        "reference": "NASA-RP-1218, Appendix D, Table D2 (published model output)",
        "reference_url": "https://ntrs.nasa.gov/citations/19890016302",
        "band_count": len(rows), "frequency_range_Hz": [FREQUENCIES[0], FREQUENCIES[-1]],
        "max_absolute_difference_dB": abs(worst["difference_dB"]),
        "max_difference_frequency_Hz": worst["frequency_Hz"],
        "rmse_dB": math.sqrt(sum(v * v for v in errors) / len(errors)),
        "mean_difference_dB": sum(errors) / len(errors),
        "peak_frequency_Hz": peak["frequency_Hz"], "peak_total_dB": peak["total_dB"],
        "regression_tolerance_dB": .01,
        "checks": ["Table D2 frequencies", "finite active sources and totals",
                   "silent disabled sources", "energy sum", "reported differences", "0.01 dB comparison"],
        "limitation": "One published-model case; not experimental validation.",
    }


def plot(rows, directory):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    frequency = [row["frequency_Hz"] for row in rows]
    fig, axes = plt.subplots(2, 1, figsize=(9, 6.6), sharex=True,
                             gridspec_kw={"height_ratios": [3, 1.3]})
    axes[0].semilogx(frequency, [r["reference_D2_dB"] for r in rows], "o", ms=5,
                     mfc="none", color="#111827", label="NASA Table D2 (model output)")
    axes[0].semilogx(frequency, [r["total_dB"] for r in rows], color="#2563eb",
                     lw=1.7, label="Aeroacoustics_cpp total")
    axes[0].set_ylabel("1/3-octave SPL (dB)")
    axes[0].set_title("BPM 1989: NACA 0012, Appendix D / Figure 45(a)\n"
                      "U = 71.3 m/s, effective angle = 1.516 deg")
    axes[0].legend(fontsize=9)
    axes[1].semilogx(frequency, [r["difference_dB"] for r in rows], ".-", color="#c2410c")
    axes[1].axhline(0, color="#6b7280", lw=.8)
    # Include zero and all errors, including when comparing a modified implementation.
    low = min(0., min(r["difference_dB"] for r in rows))
    high = max(0., max(r["difference_dB"] for r in rows))
    margin = max(.001, .1 * (high - low))
    axes[1].set_ylim(low - margin, high + margin)
    axes[1].set_ylabel("Code - D2 (dB)")
    axes[1].set_xlabel("Band center frequency (Hz)")
    for ax in axes:
        ax.grid(True, which="both", alpha=.2)
        ax.set_xlim(90, 45000)
    fig.tight_layout()
    fig.savefig(directory / "comparison.png", dpi=180)
    fig.savefig(directory / "comparison.svg")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="CSV produced by aeroacoustics_bpm1989")
    parser.add_argument("--output-dir", type=Path, default=Path("build/bpm1989-report"))
    parser.add_argument("--plot", action="store_true", help="Also plot; requires matplotlib")
    args = parser.parse_args()
    rows, summary = analyze(args.csv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    if args.plot:
        plot(rows, args.output_dir)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

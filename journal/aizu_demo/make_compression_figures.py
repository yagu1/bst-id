# journal/aizu_demo/make_compression_figures.py
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
ANALYSIS_DIR = HERE / "compression_analysis"


def read_csv(name: str):
    with (ANALYSIS_DIR / name).open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def figure_representation_size() -> None:
    rows = read_csv("compression_summary.csv")
    names = [r["representation"] for r in rows]
    values = [float(r["total_KiB"]) for r in rows]

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(names, values)
    ax.set_ylabel("Estimated total size (KiB)")
    ax.set_title("Whole-scenario representation size")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(ANALYSIS_DIR / "fig_compression_total_size.png", dpi=200)
    plt.close(fig)


def figure_record_count() -> None:
    rows = read_csv("compression_summary.csv")
    names = [r["representation"] for r in rows]
    values = [int(r["records"]) for r in rows]

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(names, values)
    ax.set_ylabel("Records")
    ax.set_yscale("log")
    ax.set_title("Whole-scenario record count")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(ANALYSIS_DIR / "fig_compression_record_count.png", dpi=200)
    plt.close(fig)


def figure_axis_merges() -> None:
    rows = read_csv("axis_merge_summary.csv")
    names = [r["axis"].upper() for r in rows]
    values = [int(r["binary_sibling_merges"]) for r in rows]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.bar(names, values)
    ax.set_xlabel("Axis")
    ax.set_ylabel("Binary sibling merges")
    ax.set_title("Normalization merges by axis")
    fig.tight_layout()
    fig.savefig(ANALYSIS_DIR / "fig_axis_merge_counts.png", dpi=200)
    plt.close(fig)


def figure_zoom_vectors() -> None:
    rows = read_csv("zoom_distribution.csv")
    rows.sort(key=lambda r: int(r["count"]), reverse=True)
    rows = rows[:20]

    labels = [
        f"({r['zx']},{r['zy']},{r['zh']},{r['zt']}) / A{r['attribute_cardinality']}"
        for r in rows
    ]
    values = [int(r["count"]) for r in rows]

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.barh(labels[::-1], values[::-1])
    ax.set_xlabel("Normalized prefix records")
    ax.set_title("Top normalized zoom vectors (X, Y, H, T)")
    fig.tight_layout()
    fig.savefig(ANALYSIS_DIR / "fig_zoom_vector_distribution.png", dpi=200)
    plt.close(fig)


def main() -> None:
    if not ANALYSIS_DIR.exists():
        raise SystemExit(
            f"{ANALYSIS_DIR} does not exist. Run analyze_compression.py first."
        )

    figure_representation_size()
    figure_record_count()
    figure_axis_merges()
    figure_zoom_vectors()

    print("Figures written to:")
    for p in sorted(ANALYSIS_DIR.glob("fig_*.png")):
        print(f"  {p}")


if __name__ == "__main__":
    main()

"""Render the single-column Figure 6 using the archived plot aggregates."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FormatStrFormatter


HERE = Path(__file__).resolve().parent
DATA = json.loads((HERE / "aggregates.json").read_text())
INK = "#203341"
MUTED = "#465867"
GRID = "#E1E7EC"
MASK = "#0072B2"
COMPLETE = "#D69500"
WIDTH_MM = 86
HEIGHT_MM = 113

plt.rcParams.update({
    "font.family": ["Arial", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "text.color": INK,
    "axes.labelcolor": INK,
    "axes.edgecolor": MUTED,
    "xtick.color": INK,
    "ytick.color": INK,
    "axes.linewidth": .5,
    "pdf.fonttype": 42,
    "axes.unicode_minus": True,
})
fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), facecolor="white")


def title(y, letter, title_text):
    fig.text(.012, y, letter, va="top", weight="bold", size=10)
    fig.text(.074, y, title_text, va="top", weight="bold", size=8.5)


def finish_axis(ax):
    ax.spines[["right", "top"]].set_visible(False)
    ax.tick_params(length=2, width=.5, pad=2)


title(.991, "A", "Response amplitude before alignment")
amplitude = fig.add_axes([.195, .775, .781, .151])
checkpoint_values = list(DATA["checkpoints"].values())
amplitude.set_axisbelow(True)
amplitude.grid(axis="y", lw=.55, color=GRID)
amplitude.vlines(range(5), 0, checkpoint_values[:5], color="#84ADB6", lw=1)
amplitude.scatter(range(5), checkpoint_values[:5], s=16, color="#306E7A", zorder=3)
amplitude.vlines(5, 0, checkpoint_values[5], color=INK, lw=1)
amplitude.scatter(5, checkpoint_values[5], s=20, marker="D", color=INK, zorder=3)
for index, value in enumerate(checkpoint_values):
    amplitude.annotate(f"{value:.4f}", (index, value), xytext=(0, 4),
                       textcoords="offset points", ha="center", va="bottom", size=7)
amplitude.axvline(4.5, color=GRID, lw=.8)
amplitude.set(xlim=(-.55, 5.55), ylim=(0, .087), ylabel="Macro SD ratio",
              xticks=range(6), xticklabels=["0", "1", "2", "3", "4", "Ensemble"],
              yticks=[0, .04, .08])
amplitude.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
finish_axis(amplitude)
fig.text(.49, .704, "Checkpoint", ha="center", size=7)

title(.662, "B", "Input alignment at fixed weights")
fig.legend(
    handles=[
        Line2D([], [], lw=3, color=MASK, label="Mask only"),
        Line2D([], [], lw=3, color=COMPLETE, label="Complete alignment"),
    ],
    loc="center left", bbox_to_anchor=(.055, .601), ncol=2, frameon=False,
    handlelength=1.1, handletextpad=.55, columnspacing=1.8, borderaxespad=0,
)

landscapes = DATA["landscape_labels"]
row_positions = list(range(len(landscapes)))
section_bottom, section_height = .204, .316
relative = fig.add_axes([.355, section_bottom, .256, section_height])
centered = fig.add_axes([.708, section_bottom, .268, section_height], sharey=relative)
for ax in (relative, centered):
    ax.set_ylim(12.64, -.64)
    ax.set_yticks([])
    ax.axhspan(11.55, 12.64, color="#EDF2F5", zorder=0)
    ax.axvline(0, color=MUTED, lw=.55)
    ax.grid(axis="x", color=GRID, lw=.5)
    ax.set_axisbelow(True)
    finish_axis(ax)
    ax.spines["left"].set_visible(False)

for name, color, offset in (
    ("mask_aligned", MASK, -.17), ("all_indices_aligned", COMPLETE, .17),
):
    condition = DATA["conditions"][name]
    y = [row + offset for row in row_positions]
    relative.barh(y, condition["sd_change_percent"], height=.28, color=color)
    centered.barh(y, [x * 1000 for x in condition["cmae_change"]],
                  height=.28, color=color)

relative.set(xlim=(-6.5, 6.5), xticks=[-6, 0, 6], xlabel="SD ratio change (%)")
centered.set(xlim=(-.23, 2.15), xticks=[0, 1, 2], xlabel="ΔcMAE (10⁻³ pKd)")
fig.text(.015, .536, "Landscape", size=7, weight="bold")
fig.text(.289, .536, "Original\nSD ratio", ha="right", size=7, weight="bold",
         linespacing=1.05)
for i, (label, baseline) in enumerate(zip(landscapes, DATA["baseline_sd_ratio"])):
    y = section_bottom + section_height * ((12.64 - i) / 13.28)
    weight = "bold" if label == "Macro" else "normal"
    fig.text(.015, y, label, va="center", weight=weight)
    ratio = f"{baseline:.4f}" if baseline >= .001 else f"{baseline:.1e}".replace("e-0", "e-")
    fig.text(.289, y, ratio, va="center", ha="right", weight=weight)

condition_columns = [.593, .767, .928]
fig.add_artist(Rectangle((.013, .012), .97, .060, transform=fig.transFigure,
                         color="#EDF2F5", linewidth=0, zorder=0))
for x, label in zip(condition_columns, ("Original", "Mask", "Complete")):
    fig.text(x, .087, label, ha="center", va="center", weight="bold")
for y, label, field in (
    (.056, "Equal residue vectors (/1055)", "equal_residue_vectors"),
    (.025, "Zero ensemble response (/211)", "output_zeros"),
):
    fig.text(.025, y, label, va="center")
    for x, condition in zip(condition_columns, DATA["conditions"].values()):
        fig.text(x, y, str(condition[field]), ha="center", va="center", weight="bold")

fig.savefig(HERE / "figure6b_source.pdf", metadata={"Creator": "Matplotlib"})

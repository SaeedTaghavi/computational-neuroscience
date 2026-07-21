"""Shared plotting style for the بیوفیزیک نورون chapters.
Matches the book aesthetic: clean white bg, light grid, English/math labels,
blue main curve, red for the 'unstable'/opposing element, green/purple trajectories.
Figures ~950 px wide (RGBA PNG).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# palette
BLUE   = "#1f77b4"
RED    = "#d62728"
GREEN  = "#2ca02c"
PURPLE = "#7e57c2"
ORANGE = "#e08214"
TEAL   = "#17a2b8"
GRAY   = "#7f7f7f"
LIGHT  = "#b0b0b0"

plt.rcParams.update({
    "figure.dpi": 130,
    "savefig.dpi": 130,
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "axes.grid": True,
    "grid.color": "#d9d9d9",
    "grid.linewidth": 0.8,
    "axes.axisbelow": True,
    "axes.edgecolor": "#444444",
    "legend.frameon": True,
    "legend.framealpha": 0.9,
    "legend.fontsize": 10.5,
    "lines.linewidth": 2.0,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.bbox": "tight",
    "savefig.transparent": False,
    "mathtext.fontset": "dejavusans",
})

def save(fig, path):
    fig.savefig(path)
    plt.close(fig)
    print("wrote", path)

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as fm

# =========================
# INPUT SETTINGS
# =========================
file_path = "Energy_breakchart.csv"

font_path = "/path/to/timesbold.ttf" #if needed

output_png = "BE.png"
output_pdf = "BE.pdf"

compounds = [
    "Control",
    "Compound1",
    "Compound2",
    "Compound3",
    "Compound4",
    "Compound5"
]

energy_term_order = [
    "van der Waal energy",
    "Electrostattic energy",
    "Polar solvation energy",
    "SASA energy",
    "Binding energy"
]

term_colors = {
    "van der Waal energy": "#1f77b4",
    "Electrostattic energy": "#2ca02c",
    "Polar solvation energy": "#8c564b",
    "SASA energy": "#7f7f7f",
    "Binding energy": "#17becf"
}

# =========================
# FONT SIZES
# =========================
AXIS_LABEL_SIZE = 24
XTICK_SIZE = 22
YTICK_SIZE = 22
LEGEND_TEXT_SIZE = 20
LEGEND_TITLE_SIZE = 20

# =========================
# FONT SETTINGS
# =========================
if os.path.isfile(font_path):
    fm.fontManager.addfont(font_path)

    axis_fp = fm.FontProperties(fname=font_path, weight="bold", size=AXIS_LABEL_SIZE)
    xtick_fp = fm.FontProperties(fname=font_path, weight="bold", size=XTICK_SIZE)
    ytick_fp = fm.FontProperties(fname=font_path, weight="bold", size=YTICK_SIZE)
    legend_fp = fm.FontProperties(fname=font_path, weight="bold", size=LEGEND_TEXT_SIZE)
    legend_title_fp = fm.FontProperties(fname=font_path, weight="bold", size=LEGEND_TITLE_SIZE)

    font_name = fm.FontProperties(fname=font_path).get_name()
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = [font_name]
else:
    axis_fp = fm.FontProperties(weight="bold", size=AXIS_LABEL_SIZE)
    xtick_fp = fm.FontProperties(weight="bold", size=XTICK_SIZE)
    ytick_fp = fm.FontProperties(weight="bold", size=YTICK_SIZE)
    legend_fp = fm.FontProperties(weight="bold", size=LEGEND_TEXT_SIZE)
    legend_title_fp = fm.FontProperties(weight="bold", size=LEGEND_TITLE_SIZE)

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# =========================
# READ CSV WITH VALUE + EB
# =========================
data_raw = pd.read_csv(file_path, header=[0, 1], sep=None, engine="python")

energy_term_col = data_raw.columns[0]
energy_terms = data_raw[energy_term_col].astype(str).str.strip()

data = data_raw.iloc[:, 1:].copy()
data.index = energy_terms

clean_columns = []
for compound, subcol in data.columns:
    clean_columns.append((str(compound).strip(), str(subcol).strip()))

data.columns = pd.MultiIndex.from_tuples(clean_columns)
data = data.apply(pd.to_numeric, errors="coerce")

available_compounds = [
    comp for comp in compounds
    if (comp, "Kcal/mol") in data.columns and (comp, "EB") in data.columns
]

available_terms = [
    term for term in energy_term_order
    if term in data.index
]

if len(available_compounds) == 0:
    raise ValueError("No valid compound columns found.")

if len(available_terms) == 0:
    raise ValueError("No valid energy terms found.")

# =========================
# PLOT
# =========================
fig, ax = plt.subplots(figsize=(14, 9))

x = np.arange(len(available_compounds))
bar_width = 0.15

for i, term in enumerate(available_terms):
    means = []
    ebs = []

    for comp in available_compounds:
        means.append(data.loc[term, (comp, "Kcal/mol")])
        ebs.append(data.loc[term, (comp, "EB")])

    ax.bar(
        x + i * bar_width,
        means,
        width=bar_width,
        label=term,
        color=term_colors.get(term, None),
        yerr=ebs,
        capsize=5,
        edgecolor="none",
        linewidth=0,
        error_kw={
            "elinewidth": 1.8,
            "capthick": 1.8,
            "ecolor": "black"
        }
    )

# Zero line
ax.axhline(0, color="none", linewidth=1.4)

# Axis labels
ax.set_xlabel("Compounds", fontproperties=axis_fp, labelpad=15)
ax.set_ylabel("Energy Values (Kcal/mol)", fontproperties=axis_fp, labelpad=15)

# X tick labels
ax.set_xticks(x + bar_width * (len(available_terms) - 1) / 2)
ax.set_xticklabels(
    available_compounds,
    rotation=0,
    ha="right"
)

# Force x tick font
for label in ax.get_xticklabels():
    label.set_fontproperties(xtick_fp)

# Force y tick font
for label in ax.get_yticklabels():
    label.set_fontproperties(ytick_fp)

# Tick line size
ax.tick_params(
    axis="both",
    which="major",
    direction="out",
    length=8,
    width=2.2
)

# Legend in top-middle
legend = ax.legend(
    title="Energy Terms",
    loc="upper center",
    bbox_to_anchor=(0.45, 1.02),
    frameon=False,
    ncol=1,
    prop=legend_fp,
    title_fontproperties=legend_title_fp,
    handlelength=2.8,
    handleheight=1.0,
    labelspacing=0.4
)

# Box line width
for spine in ["bottom", "left", "top", "right"]:
    ax.spines[spine].set_linewidth(2.2)

# Axis limits
ax.set_ylim(-65, 48)

plt.tight_layout()

plt.savefig(output_png, dpi=600, bbox_inches="tight")
plt.savefig(output_pdf, dpi=600, bbox_inches="tight")

plt.show()

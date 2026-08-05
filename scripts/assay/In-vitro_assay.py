# cetp_inhibition_torc100_no_scaling.py
# Plots % inhibition bars: Torcetrapib (fixed 100%), then Adapalene 1/10, Buclizine 1/10, Mefloquine 1/10.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

INFILE = "Inhibitor_assay.csv"   # your table with 3 replicate "% Inhibition" columns
OUT_FIG = "Inhibitor_assay.png"
OUT_SUM = "Inhibitor_assay.csv"

# --- robust read (handles µ in headers / Windows encodings) ---
def read_csv_robust(path):
    for enc in ("utf-8", "utf-8-sig", "cp1252", "latin-1"):
        try:
            return pd.read_csv(path, encoding=enc, engine="python")
        except UnicodeDecodeError:
            continue
    return pd.read_csv(path, encoding="latin-1", engine="python", encoding_errors="replace")

df = read_csv_robust(INFILE)
df.columns = [c.strip().replace("µ","u").replace("μ","u") for c in df.columns]

# column ids
cmp_col  = "Compound_name" if "Compound_name" in df.columns else "Compound"
dose_col = "Concentration (uM)" if "Concentration (uM)" in df.columns else (
           "Concentration_uM" if "Concentration_uM" in df.columns else "Concentration")

# replicate columns (must start with '% Inhibition')
rep_cols = [c for c in df.columns if str(c).strip().lower().startswith("% inhibition")]
if not rep_cols:
    raise ValueError("Input must contain three columns starting with '% Inhibition'.")

# numeric
df[dose_col] = pd.to_numeric(df[dose_col], errors="coerce")
for c in rep_cols: df[c] = pd.to_numeric(df[c], errors="coerce")

# long → summary (mean ± SD)
long = df.melt(id_vars=[cmp_col, dose_col], value_vars=rep_cols,
               var_name="Rep", value_name="%Inh").dropna()
summ = (long.groupby([cmp_col, dose_col], as_index=False)["%Inh"]
        .agg(Mean="mean", SD="std"))

# helper to pull a bar
def bar(drug, dose):
    row = summ[(summ[cmp_col].str.strip()==drug) & np.isclose(summ[dose_col], dose)]
    if row.empty: return np.nan, np.nan
    return float(row["Mean"].values[0]), float(row["SD"].values[0])

# assemble bars (torcetrapib fixed to 100% without scaling others)
labels = ["Control",
          "Drug 1 conc1","Drug 1 conc2",
          "Drug 2 conc1","Drug 2 conc2",
          "Drug 3 conc1","Drug 3 conc2"]

means  = [100.0]  # torcetrapib fixed 100%
errs   = [0.0]

for drug, dose in [("Drug 1",conc1), ("Drug 1",conc2),
                   ("Drug 2",conc1), ("Drug 2",conc2),
                   ("Drug 3",conc1), ("Drug 3",conc2)]:
    m, s = bar(drug, dose)
    means.append(m)
    errs.append(s)

# save the exact numbers used
pd.DataFrame({"Label": labels, "%Inhibition": means, "SD": errs}).to_csv(OUT_SUM, index=False)

# colors: control black; blues for Adapalene; greens for Buclizine & Mefloquine
COL_TORC = "#000000"
BLUES    = ["#77AFC4", "#5A7798"]
GREENS_B = ["#9CC6B1", "#76A78F"]
GREENS_M = ["#C6DEC9", "#93C0A1"]
colors = [COL_TORC] + BLUES + GREENS_B + GREENS_M

# plot
fig, ax = plt.subplots(figsize=(8, 7))
x = np.arange(len(labels))
ax.bar(x, means, yerr=errs, capsize=3, color=colors)
# Increase spine (box) line width and push ticks outside
for spine in ["bottom", "left"]:
    ax.spines[spine].set_linewidth(2.0)   # thick axis box lines
for spine in ["top", "right"]:
    ax.spines[spine].set_linewidth(2.0)   # hide opposite box lines
ax.set_ylabel("% inhibition of protein activity", fontsize=16)
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=14, rotation=55, ha="right")
ax.tick_params(axis='y', labelsize=14)   # << increase y-axis tick label size
ax.set_ylim(0, max(110, np.nanmax(means)*1.15))
ax.set_title("Protein inhibition: Control + test compounds")
ax.grid(False)

plt.tight_layout()
plt.savefig(OUT_FIG, dpi=600)
plt.show()
print(f"Saved: {OUT_FIG} and {OUT_SUM}")

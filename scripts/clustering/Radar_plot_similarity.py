from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import font_manager
from matplotlib.colors import TwoSlopeNorm, Normalize
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

# ========================= USER SETTINGS =========================
CSV = "table.csv"      # your similarity table
ID_COL = "Ligand_name"            # ligand identifier column
K = 30                            # number of clusters
N_RADAR = 6                       # show top 6 largest clusters
RANDOM_STATE = 7
DPI = 600

OUTDIR = Path("Figure3_outputs")
OUTDIR.mkdir(exist_ok=True, parents=True)

FONT_PATH = "timesbold.ttf"       # set None if not needed

USE_ZSCORE = True                 # recommended
SHOW_COLORBAR = True              # set False if you want no colorbar
# ================================================================

# -------- metric columns --------
tmp = pd.read_csv(CSV)
cols = set(tmp.columns)

if "Russell_similarity" in cols:
    russell_col = "Russell_similarity"
elif "Russel_similarity" in cols:
    russell_col = "Russel_similarity"
else:
    raise ValueError("Could not find Russell/Russel similarity column.")

TWO_D_METRICS = [
    "Tanimoto_similarity",
    "Cosine_similarity",
    "Dice_similarity",
    "Kulczynski_similarity",
    russell_col,
    "Sokal_similarity",
    "Tversky_similarity",
]
THREE_D_METRICS = ["USR_similarity"]
ELEC_METRICS = ["MEP_similarity"]

METRICS = TWO_D_METRICS + THREE_D_METRICS + ELEC_METRICS

GROUPS = {
    "2D fingerprint": TWO_D_METRICS,
    "3D shape": THREE_D_METRICS,
    "Electrostatic": ELEC_METRICS
}

GROUP_COLORS = {
    "2D fingerprint": "#66C2A5",   # teal-green
    "3D shape": "#8DA0CB",         # soft blue-purple
    "Electrostatic": "#FC8D62"     # soft orange
}

POLY_COLOR = "#2F7F6F"
BG_ALPHA = 0.10

# -------- font --------
if FONT_PATH and Path(FONT_PATH).is_file():
    font_manager.fontManager.addfont(FONT_PATH)
    fam = font_manager.FontProperties(fname=FONT_PATH).get_name()
    plt.rcParams["font.family"] = fam

plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "font.weight": "bold",
    "axes.linewidth": 1.0,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.transparent": False,
})

# -------- load data --------
df = pd.read_csv(CSV)
required = [ID_COL] + METRICS
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

df = df.dropna(subset=required).copy()
X = df[METRICS].astype(float).to_numpy()

# -------- standardize --------
if USE_ZSCORE:
    scaler = StandardScaler()
    M = scaler.fit_transform(X)
    cbar_label = "Cluster centroid (z-score)"
else:
    M = X.copy()
    cbar_label = "Descriptor value"

# -------- clustering --------
km = KMeans(n_clusters=K, n_init=30, random_state=RANDOM_STATE)
df["cluster_k30"] = km.fit_predict(M) + 1

cluster_sizes = df["cluster_k30"].value_counts().sort_values(ascending=False)
top_clusters = cluster_sizes.head(N_RADAR).index.tolist()

# -------- centroids --------
cluster_centroids = {}
all_centroid_values = []

for clus in top_clusters:
    sub = M[df["cluster_k30"].values == clus]
    centroid = sub.mean(axis=0)
    cluster_centroids[clus] = centroid
    all_centroid_values.extend(centroid.tolist())

all_centroid_values = np.array(all_centroid_values)

# -------- color normalization --------
vmax = np.max(np.abs(all_centroid_values))
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
cmap = plt.cm.RdYlBu_r

# -------- helper: sector background --------
def add_sector_background(ax, theta, rmax):
    """
    Add light background sectors for the 3 descriptor classes.
    """
    sector_defs = [
        ("2D fingerprint", 0, len(TWO_D_METRICS)-1),
        ("3D shape", len(TWO_D_METRICS), len(TWO_D_METRICS)),
        ("Electrostatic", len(TWO_D_METRICS)+1, len(TWO_D_METRICS)+1),
    ]

    delta = theta[1] - theta[0]
    for name, i0, i1 in sector_defs:
        start = theta[i0] - delta/2
        end = theta[i1] + delta/2
        tt = np.linspace(start, end, 200)

        # light filled sector
        ax.fill_between(
            tt,
            0,
            rmax,
            color=GROUP_COLORS[name],
            alpha=BG_ALPHA,
            zorder=0
        )

        # outer arc
        ax.plot(tt, np.full_like(tt, rmax*0.995),
                color=GROUP_COLORS[name], lw=3.0, solid_capstyle="round", zorder=3)

        # sector label
        mid = 0.5 * (start + end)
        ax.text(mid, rmax*1.11, name,
                ha="center", va="center",
                fontsize=10, fontweight="bold",
                color=GROUP_COLORS[name])

# -------- helper: radar --------
def draw_radar(ax, values, cluster_id, n_members):
    L = len(values)
    theta = np.linspace(0, 2*np.pi, L, endpoint=False)

    vals_closed = np.r_[values, values[0]]
    theta_closed = np.r_[theta, theta[0]]

    # radial limits
    vmin = np.min(all_centroid_values)
    vmax = np.max(all_centroid_values)
    span = (vmax - vmin) if vmax > vmin else 1.0
    pad = 0.20 * span
    rmin = min(vmin - pad, -0.1)
    rmax = vmax + pad

    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_ylim(rmin, rmax)

    # background sectors
    add_sector_background(ax, theta, rmax)

    # polygon
    ax.plot(theta_closed, vals_closed, color=POLY_COLOR, lw=2.0, zorder=4)
    ax.fill(theta_closed, vals_closed, color=POLY_COLOR, alpha=0.14, zorder=2)

    # colored points
    ax.scatter(theta, values, c=values, cmap=cmap, norm=norm,
               s=38, edgecolors="#222222", linewidths=0.35, zorder=5)

    # no individual spoke labels
    ax.set_xticks(theta)
    ax.set_xticklabels([""] * L)

    # hide radial tick labels
    ax.set_yticks([])

    # spoke grid
    ax.grid(color="#B5B5B5", alpha=0.75, linewidth=0.8)

    # title
    ax.set_title(f"Cluster {cluster_id} (n={n_members})",
                 y=1.12, fontsize=12, fontweight="bold")

# -------- plot layout --------
ncols = 3
nrows = int(np.ceil(N_RADAR / ncols))

fig, axes = plt.subplots(
    nrows, ncols,
    subplot_kw=dict(polar=True),
    figsize=(13.5, 8.0),
    constrained_layout=True
)
axes = np.array(axes).reshape(-1)

for ax, clus in zip(axes, top_clusters):
    centroid = cluster_centroids[clus]
    n_members = int(cluster_sizes[clus])
    draw_radar(ax, centroid, clus, n_members)

for ax in axes[len(top_clusters):]:
    ax.axis("off")

# color bar
if SHOW_COLORBAR:
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.tolist(), fraction=0.024, pad=0.04)
    cbar.set_label(cbar_label, fontsize=12)

fig.suptitle("Cluster profiles across descriptor classes",
             fontsize=16, fontweight="bold", y=1.02)

out_png = OUTDIR / "Radar_plot.png"
out_pdf = OUTDIR / "Radar_plot.pdf"

fig.savefig(out_png, dpi=DPI, bbox_inches="tight")
fig.savefig(out_pdf, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")

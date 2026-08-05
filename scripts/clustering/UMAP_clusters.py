from pathlib import Path
import os
os.environ["MPLBACKEND"] = "Agg"          # non-GUI backend for Matplotlib
os.environ["QT_QPA_PLATFORM"] = "offscreen"

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")                            # belt-and-suspenders
import matplotlib.pyplot as plt
from matplotlib import font_manager
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score
from scipy.spatial import ConvexHull
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# ---------------- user knobs ----------------
CSV = "File"   # has columns: Ligand_name + similarity metrics
ID  = "Ligand_name"
METRICS = ["Tanimoto_similarity","Cosine_similarity","Dice_similarity",
           "Kulczynski_similarity","Russel_similarity","Sokal_similarity",
           "Tversky_similarity","USR_similarity","MEP_similarity"]

K = 30
METHOD = "umap"        # "umap" (preferred) or "tsne"
UMAP_NEIGHBORS = 30
UMAP_MIN_DIST  = 0.12
TSNE_PERPLEXITY = 40
TSNE_LR         = 240
TSNE_EXAGG      = 1.2

DRAW_HULLS   = True
MIN_HULL_SIZE = 6      # don't draw a hull for very tiny clusters
POINT_SIZE   = 18
POINT_ALPHA  = 0.85
OUT          = "UMAP_cluster.png"
DPI          = 600
# --------------------------------------------

# Fonts: Times Bold
font_manager.fontManager.addfont("timesbold.ttf")
plt.rcParams["font.family"] = font_manager.FontProperties(fname="timesbold.ttf").get_name()
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.linewidth"] = 1.1

# Load & z-score
df = pd.read_csv(CSV)
Z  = StandardScaler().fit_transform(df[METRICS].values)

# If you already have labels, use them; else fit KMeans here
if "cluster_k30_min10" in df.columns:
    labels = df["cluster_k30_min10"].to_numpy()
elif "cluster_k30" in df.columns:
    labels = df["cluster_k30"].to_numpy()
else:
    km = KMeans(n_clusters=K, n_init=30, random_state=7).fit(Z)
    labels = km.labels_ + 1
    df["cluster_k30"] = labels

# Embedding
if METHOD.lower() == "umap":
    try:
        import umap
        reducer = umap.UMAP(n_neighbors=UMAP_NEIGHBORS, min_dist=UMAP_MIN_DIST,
                            metric="euclidean", random_state=7)
        XY = reducer.fit_transform(Z)
        title = f"UMAP of 9-metric space colored by clusters (k={len(np.unique(labels))})"
    except Exception as e:
        print("UMAP not available, falling back to t-SNE. Error:", e)
        METHOD = "tsne"
if METHOD.lower() == "tsne":
    tsne = TSNE(n_components=2, perplexity=TSNE_PERPLEXITY,
                learning_rate=TSNE_LR, init="pca",
                early_exaggeration=TSNE_EXAGG, random_state=7, n_iter=1200)
    XY = tsne.fit_transform(Z)
    title = f"t-SNE of 9-metric space colored by clusters (k={len(np.unique(labels))})"

# Palette (30 distinct colors) – use non-deprecated API
palette = (list(mpl.colormaps["tab20"].colors) + list(mpl.colormaps["Set3"].colors))[:len(np.unique(labels))]
color_map = {k: palette[i % len(palette)] for i, k in enumerate(sorted(np.unique(labels)))}

# Plot
fig, ax = plt.subplots(figsize=(8.2, 6.8), constrained_layout=True)

# Optional convex hulls (light)
if DRAW_HULLS:
    for k in sorted(np.unique(labels)):
        idx = np.where(labels == k)[0]
        if len(idx) >= MIN_HULL_SIZE:
            pts = XY[idx, :]
            try:
                hull = ConvexHull(pts)
                poly = pts[hull.vertices]
                ax.fill(poly[:,0], poly[:,1],
                        color=color_map[k], alpha=0.12, linewidth=0, zorder=1)
                # thin outline helps readability
                ax.plot(np.r_[poly[:,0], poly[0,0]], np.r_[poly[:,1], poly[0,1]],
                        color=color_map[k], lw=0.5, alpha=0.35, zorder=2)
            except Exception:
                pass

# Points
for k in sorted(np.unique(labels)):
    idx = np.where(labels == k)[0]
    ax.scatter(XY[idx,0], XY[idx,1], s=POINT_SIZE,
               color=color_map[k], alpha=POINT_ALPHA,
               edgecolors="white", linewidths=0.35, rasterized=True, zorder=3)

# Label cluster centroids with numbers
for k in sorted(np.unique(labels)):
    idx = np.where(labels == k)[0]
    cx, cy = XY[idx,0].mean(), XY[idx,1].mean()
    ax.text(cx, cy, str(k), ha="center", va="center",
            fontsize=10, color="black",
            bbox=dict(boxstyle="circle,pad=0.28", fc="white", ec="black", lw=0.8, alpha=0.95),
            zorder=4)

ax.set_xlabel("Component 1", fontsize=14)
ax.set_ylabel("Component 2", fontsize=14)
ax.tick_params(axis="both", labelsize=14)
ax.set_title(title, fontsize=14)
ax.grid(False)
for s in ["top","right"]:
    ax.spines[s].set_visible(False)

# ---- Silhouette inset (cluster-size bars + value) ----
try:
    sil = silhouette_score(Z, labels)
except Exception:
    sil = float("nan")
sizes = pd.Series(labels).value_counts().sort_values(ascending=False).values
ax_in = inset_axes(ax, width="22%", height="18%", loc="upper right", borderpad=1.0)
ax_in.bar(range(len(sizes)), sizes, color="#bdbdbd", width=0.9)
ax_in.set_title(f"Silhouette={sil:.3f}" if np.isfinite(sil) else "Silhouette=N/A", fontsize=12)
ax_in.set_xticks([]); ax_in.set_yticks([])
for sp in ax_in.spines.values(): sp.set_visible(False)

plt.savefig(OUT, dpi=DPI, bbox_inches="tight")
print("Saved", OUT)

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib import cm
from scipy.cluster import hierarchy as sch
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.manifold import TSNE

# ------------------- user knobs -------------------
CSV = "file.csv"     # your file
ID  = "Ligand_name"
METRICS = [
    "Tanimoto_similarity", "Cosine_similarity", "Dice_similarity",
    "Kulczynski_similarity", "Russel_similarity", "Sokal_similarity",
    "Tversky_similarity", "USR_similarity", "MEP_similarity"
]
OUTDIR = Path("sim_figs2"); OUTDIR.mkdir(exist_ok=True)
K = 30                    
N_RADAR = 6              
TOP_PARALLEL = 30        
RANDOM_STATE = 7

# fonts
font_manager.fontManager.addfont("timesbold.ttf")
plt.rcParams["font.family"] = font_manager.FontProperties(fname="timesbold.ttf").get_name()
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.linewidth"] = 1.1

# ------------------- load & prep -------------------
df = pd.read_csv(CSV)
X = df[METRICS].astype(float).values
ids = df[ID].astype(str).values

# standardize metrics (z-score)
scaler = StandardScaler()
Z = scaler.fit_transform(X)

# consensus score: mean of z across metrics (higher = more similar by consensus)
consensus = Z.mean(axis=1)
df["consensus_zmean"] = consensus
df["consensus_rank"] = (-consensus).argsort().argsort() + 1  # 1 = best

# -------------- Figure A: metric agreement --------------
# Spearman correlations between metric columns across compounds
R = np.zeros((len(METRICS), len(METRICS)))
for i, a in enumerate(METRICS):
    for j, b in enumerate(METRICS):
        R[i, j] = pd.Series(df[a]).corr(df[b], method="spearman")

# cluster metrics for a pleasing order
link = sch.linkage(1 - R, method="average")  # 1-R ~ dissimilarity
order = sch.leaves_list(link)
Rr = R[order][:, order]
labels_ordered = [METRICS[i] for i in order]

fig, ax = plt.subplots(figsize=(6.6, 6.0), constrained_layout=True)
im = ax.imshow(Rr, vmin=0, vmax=1, cmap="Greys")
ax.set_xticks(range(len(METRICS)), labels_ordered, rotation=45, ha="right", fontsize=14)
ax.set_yticks(range(len(METRICS)), labels_ordered, fontsize=14)
ax.set_title("Agreement between similarity metrics (Spearman ρ)", fontsize=14)
cbar = fig.colorbar(im, ax=ax); cbar.set_label("ρ")
fig.savefig(OUTDIR/"FigA_metric_agreement.png", dpi=600, bbox_inches="tight")
plt.close(fig)

# -------------- Figure B: PCA biplot --------------
pca = PCA(n_components=2, random_state=RANDOM_STATE)
XY = pca.fit_transform(Z)
expl = pca.explained_variance_ratio_
fig, ax = plt.subplots(figsize=(7.2, 6.0), constrained_layout=True)
ax.scatter(XY[:,0], XY[:,1], s=12, c="#8DA0CB", alpha=0.65, edgecolors="none", rasterized=True)
ax.set_xlabel(f"PC1 ({expl[0]*100:.1f}%)", fontsize=14)
ax.set_ylabel(f"PC2 ({expl[1]*100:.1f}%)", fontsize=14)
ax.set_title("PCA of similarity metrics (biplot)", fontsize=14)

# metric loadings as arrows
load = pca.components_.T  # 9×2
scale = 2.3
for i, m in enumerate(METRICS):
    ax.arrow(0, 0, load[i,0]*scale, load[i,1]*scale, color="black", head_width=0.07, length_includes_head=True)
    ax.text(load[i,0]*scale*1.1, load[i,1]*scale*1.1, m.replace("_similarity",""), fontsize=12)

ax.axhline(0, color="black", lw=0.8); ax.axvline(0, color="black", lw=0.8)
fig.savefig(OUTDIR/"FigB_PCA_biplot.png", dpi=600, bbox_inches="tight")
plt.close(fig)

# -------------- Figure C: consensus score + rank concordance --------------
fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)

# histogram of consensus
axes[0].hist(consensus, bins=40, color="#68c3ad", alpha=0.9)
axes[0].set_title("Consensus (z-mean across 9 metrics)", fontsize=14)
axes[0].set_xlabel("Consensus score (higher = better)", fontsize=14)
axes[0].set_ylabel("Count", fontsize=14)

# Spearman ρ between each metric and consensus
rho = [pd.Series(df[m]).corr(df["consensus_zmean"], method="spearman") for m in METRICS]
axes[1].barh(range(len(METRICS)), rho, color="#bdbdbd")
axes[1].set_yticks(range(len(METRICS)), [m.replace("_similarity","") for m in METRICS])
axes[1].set_xlabel("Spearman ρ with consensus")
axes[1].set_xlim(0, 1)
axes[1].invert_yaxis()
axes[1].set_title("Metric–consensus agreement")

fig.savefig(OUTDIR/"FigC_consensus_and_agreement.png", dpi=600, bbox_inches="tight")
plt.close(fig)

# -------------- Figure D: K-means, elbow + silhouette, then embed & sizes --------------
# quick elbow & silhouette for k grid
k_list = list(range(10, 91, 10))
inertias, sils = [], []
for k in k_list:
    km = KMeans(n_clusters=k, n_init=20, random_state=RANDOM_STATE)
    labs = km.fit_predict(Z)
    inertias.append(km.inertia_)
    # silhouette on Euclidean in z-space
    sils.append(silhouette_score(Z, labs))

fig, ax = plt.subplots(figsize=(6.0, 4.2), constrained_layout=True)
ax.plot(k_list, inertias, "-o")
ax.set_xlabel("k", fontsize=14); ax.set_ylabel("K-means inertia", fontsize=14); ax.set_title("Elbow scan", fontsize=14)
fig.savefig(OUTDIR/"FigD1_elbow.png", dpi=600, bbox_inches="tight"); plt.close(fig)

fig, ax = plt.subplots(figsize=(6.0, 4.2), constrained_layout=True)
ax.plot(k_list, sils, "-o")
ax.axvline(K, color="black", ls="--", lw=1.2)
ax.set_xlabel("k", fontsize=14); ax.set_ylabel("Silhouette", fontsize=14); ax.set_title("Silhouette vs k", fontsize=14)
fig.savefig(OUTDIR/"FigD2_silhouette.png", dpi=600, bbox_inches="tight"); plt.close(fig)

# final k-means with k=40
km = KMeans(n_clusters=K, n_init=30, random_state=RANDOM_STATE)
labels = km.fit_predict(Z)
df["cluster_k30"] = labels + 1  # 1..K

# 2-D t-SNE on z-space (fast & clear)
tsne = TSNE(n_components=2, perplexity=30, learning_rate=200,
            n_iter=1000, init="random", random_state=RANDOM_STATE)
xy = tsne.fit_transform(Z)

# color palette
colors = cm.get_cmap("tab20").colors
fig, ax = plt.subplots(figsize=(7.2, 6.0), constrained_layout=True)
for k in range(1, K+1):
    mask = df["cluster_k30"] == k
    ax.scatter(xy[mask,0], xy[mask,1], s=10, color=colors[(k-1) % len(colors)],
               alpha=0.8, edgecolors="none", rasterized=True)
ax.set_xlabel("t-SNE 1", fontsize=14); ax.set_ylabel("t-SNE 2", fontsize=14)
ax.set_title(f"t-SNE of 9-metric space colored by K-means (k={K})", fontsize=14)
fig.savefig(OUTDIR/"FigD3_tsne_k30.png", dpi=600, bbox_inches="tight"); plt.close(fig)

# cluster sizes
sizes = df["cluster_k30"].value_counts().sort_index()
fig, ax = plt.subplots(figsize=(9.0, 4.0), constrained_layout=True)
ax.bar(range(1, K+1), sizes.values, color="#8DA0CB")
ax.set_xlabel("Cluster #", fontsize=14); ax.set_ylabel("Size", fontsize=14)
ax.set_title("Cluster size distribution (k=30)", fontsize=14)
fig.savefig(OUTDIR/"FigD4_cluster_sizes.png", dpi=600, bbox_inches="tight"); plt.close(fig)

# -------------- Figure E: radar plots for a few cluster centroids --------------
def radar(ax, values, labels, color="#2f7f6f"):
    # closed polygon
    L = len(values)
    theta = np.linspace(0, 2*np.pi, L, endpoint=False)
    values = np.r_[values, values[0]]
    theta  = np.r_[theta, theta[0]]
    ax.plot(theta, values, color=color, lw=2)
    ax.fill(theta, values, color=color, alpha=0.20)
    ax.set_xticks(np.linspace(0, 2*np.pi, L, endpoint=False), labels, fontsize=12)
    ax.set_yticks([])
    ax.set_ylim(-2.0, 2.0)  # z-score range

# pick N_RADAR largest clusters
top_k = sizes.sort_values(ascending=False).head(N_RADAR).index.tolist()
fig, axes = plt.subplots(2, (N_RADAR+1)//2, subplot_kw=dict(polar=True),
                         figsize=(11, 6.0), constrained_layout=True)
axes = np.array(axes).reshape(-1)
labels_clean = [m.replace("_similarity","") for m in METRICS]
for ax, k in zip(axes, top_k):
    sub = Z[df["cluster_k30"] == k]  # already z-scored
    centroid = sub.mean(axis=0)
    radar(ax, centroid, labels_clean)
    ax.set_title(f"Cluster {k} (n={len(sub)})", fontsize=14, y=1.08)
# hide unused axes
for ax in axes[len(top_k):]:
    ax.axis("off")
fig.suptitle("Cluster “fingerprints” (z-scored metric centroids)", fontsize=14)
fig.savefig(OUTDIR/"FigE_radar_centroids.png", dpi=600, bbox_inches="tight"); plt.close(fig)

# -------------- Figure F: parallel coordinates for top compounds --------------
top = df.sort_values("consensus_zmean", ascending=False).head(TOP_PARALLEL)
Y = (top[METRICS] - top[METRICS].mean()) / top[METRICS].std(ddof=0)  # z within full set
fig, ax = plt.subplots(figsize=(10.5, 5.2), constrained_layout=True)
xpos = np.arange(len(METRICS))
for i, row in Y.iterrows():
    ax.plot(xpos, row.values, color=(0,0,0,0.15), lw=1.2)
ax.set_xticks(xpos, labels_clean, rotation=0, fontsize=14)
ax.set_ylabel("z-score per metric", fontsize=14)
ax.set_title(f"Parallel coordinates — top {TOP_PARALLEL} by consensus", fontsize=14)
ax.axhline(0, color="black", lw=0.8)
fig.savefig(OUTDIR/"FigF_parallel_top.png", dpi=600, bbox_inches="tight"); plt.close(fig)

print("Saved figures in:", OUTDIR.resolve())

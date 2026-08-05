import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# 1) Register your font and set it *once*, *before* any plotting or rcParams tweaks:
font_path = "/path/to/timesbold.ttf"
fm.fontManager.addfont(font_path)
prop = fm.FontProperties(fname=font_path)
font_name = prop.get_name()

plt.rcParams['font.family']      = 'serif'
plt.rcParams['font.serif']       = [font_name]
plt.rcParams['font.weight']      = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
import seaborn as sns

# Tell seaborn not to rescale your font sizes
sns.set_context("notebook", rc={
    "axes.titlesize": 18,
    "axes.labelsize": 16,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18
})
sns.set_style("white")


# 2) Now do your Seaborn styling and plot:
sns.set(style="white")

df = pd.read_csv("BindingEnergy_vs.csv")
fig, ax = plt.subplots(figsize=(10, 6))

sns.violinplot(
    x="Compounds",
    y="Binding Energy (Kcal/mol)",
    data=df,
    palette="muted",
    ax=ax
)

# 3) Set titles and labels *after* all rcParams are in place:
#ax.set_title("Binding Energy Distribution by Compound",
#             fontproperties=prop)
ax.set_xlabel("Compounds",
              fontproperties=prop)
ax.set_ylabel("Binding Energy (Kcal/mol)",
              fontproperties=prop)
#ax.title.set_size(18)
ax.xaxis.label.set_size(18)
ax.yaxis.label.set_size(18)

# 4) Ticks – enforce both your font and size:
for tl in ax.get_xticklabels() + ax.get_yticklabels():
    tl.set_fontproperties(prop)
    tl.set_fontsize(16)           # explicit tick-label size
    tl.set_fontweight('bold')

# 5) Border thickness:
for spine in ax.spines.values():
    spine.set_linewidth(2)

plt.tight_layout()
plt.savefig("violin_plot_custom.png", dpi=600, bbox_inches="tight")
plt.show()

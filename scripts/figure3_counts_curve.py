import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

args = sys.argv[1:]
print(args)

palette = {
    "MERFISH": "#0066CC",
    "EELFISH": "#FC3B11",
    "Resolve": "#DE1222",
    "10x": "#59C230",
    "Vizgen": "#E96A9E",
}
cutoffs = [0, 0.25, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.97, 0.98, 0.99]
ndatasets = 5

# molecules
mols = [pd.read_csv(arg) for arg in args[:ndatasets]]
mols = [mol[pd.notna(mol['cell'])] for mol in mols]

ncell = [len(mol["cell"].unique()) for mol in mols]

# annotations
anns = [pd.read_csv(arg) for arg in args[ndatasets:(ndatasets * 2)]]

# load in and filter markers
markers = pd.read_csv(args[ndatasets * 2])
if (args[-3] == "True"):
    markers = markers[(markers["avg_log2FC"] > 0.5) & (markers["pct.2"] < 0.05)]
else:
    markers = markers[(markers["avg_log2FC"] > 0.5) & (markers["pct.2"] < 0.25)]

genes = [mol.gene.unique() for mol in mols]
shared_genes = genes[0]
for gene in genes[1:]:
    shared_genes = np.intersect1d(shared_genes, gene)

marker_genes = markers.gene.unique()
shared_marker_genes = np.intersect1d(shared_genes, marker_genes)

celltype_marker = {}
celltypes_covered = []
cts = anns[0]["predicted.class"].unique()

for ann in anns[1:]:
    cts = np.intersect1d(cts, ann["predicted.class"].unique())

if (args[-3] == "True"):
    for g in marker_genes:
        ss = markers[markers.gene == g]
        m = max(ss.avg_log2FC)
        ct = ss[ss.avg_log2FC == m].cluster.iloc[0]
        if ct in cts:
            celltype_marker[g] = ct
            celltypes_covered.append(ct)
    celltypes_covered = set(celltypes_covered)
else:
    for g in shared_marker_genes:
        ss = markers[markers.gene == g]
        m = max(ss.avg_log2FC)
        ct = ss[ss.avg_log2FC == m].cluster.iloc[0]
        if ct in cts:
            celltype_marker[g] = ct
            celltypes_covered.append(ct)
    celltypes_covered = set(celltypes_covered)

# construct gene expression matrices for each dataset at various baysor cutoffs
gex_matrices = [[] for _ in range(ndatasets)]
for cutoff in cutoffs:
    for i in range(ndatasets):
        df_ss = mols[i][mols[i]['assignment_confidence'] >= cutoff]
        ct = pd.crosstab(index=df_ss["gene"], columns=df_ss["cell"])
        gex_matrices[i].append(ct)

for i in range(ndatasets):
    if (args[-3] == "True"):
        mols[i] = mols[i][mols[i].gene.isin(marker_genes)]
    else:
        mols[i] = mols[i][mols[i].gene.isin(shared_marker_genes)]
    d = {k: v for k, v in zip(anns[i]["cell"], anns[i]["predicted.class"])}
    mols[i]['class'] = mols[i]['cell'].map(d)
    mols[i] = mols[i][mols[i]["class"].isin(celltypes_covered)]
    mols[i]['gene_class'] = mols[i]['gene'].map(celltype_marker)

ncells = [len(mol["cell"].unique()) for mol in mols]

matrices = [[] for _ in range(ndatasets)]

for cutoff in cutoffs:
    for i in range(ndatasets):
        df_ss = mols[i][mols[i]['assignment_confidence'] >= cutoff]
        ct = pd.crosstab(index=df_ss["class"], columns=df_ss["gene_class"])
        matrices[i].append(ct)

print(f"SHARED MARKER GENES: {shared_marker_genes}")

mecr_ref_markers = pd.read_csv(args[ndatasets * 2], index_col=0)
ref_markers_ss = mecr_ref_markers[(mecr_ref_markers['pct.1'] > 0.5) & (mecr_ref_markers['pct.2'] < 0.1)]
ref_markers_ss = ref_markers_ss.drop_duplicates(subset='gene', keep=False)

def get_coexpression_rate(obj, sc_markers, genes):
    coexp_rates = []
    mtx = obj.loc[genes]
    for g1 in genes:
        for g2 in genes:
            if g1 != g2 and g1 > g2 and sc_markers.at[g1, 'cluster'] != sc_markers.at[g2, 'cluster']:
                c1 = mtx.loc[g1]
                c2 = mtx.loc[g2]
                rate = np.sum((c1 > 0) & (c2 > 0)) / np.sum((c1 > 0) | (c2 > 0))
                coexp_rates.append(rate)
    return coexp_rates

### Generate TP and FP data frame for each cutoff
tp, fp, tech, color, ncells, mean_fcs, med_cs, med_fcs, threshold_mean_cs = [], [], [], [], [], [], [], [], []
for i in range(ndatasets):
    for m in matrices[i]:
        tech.append(args[(ndatasets * 2) + i + 1])
        # add rep2 to palette
        idx = args[(ndatasets * 2) + i + 1].split("_")[0]
        palette[args[(ndatasets * 2) + i + 1]] = palette[idx]
        color.append(palette[idx])
    # compute MECR
    for gm in gex_matrices[i]:
        j = len(cutoffs) - 1
        if (len(tech) <= 2 and tech[-1] != tech[-2]):
            j += len(cutoffs)
        mecr = np.mean(
            get_coexpression_rate(gm, ref_markers_ss, list(set(gex_matrices[i][-1].index).intersection(set(ref_markers_ss.index))))
        )
        fp.append(mecr)
        # compute the average gene counts per cell in the gm matrix and append to tp
        print(gm.shape)
        print(tech[-1])
        cell_counts = np.sum(gm, axis=0)
        feature_counts = np.sum(gm > 0, axis=0)
        ncells.append(gm.shape[-1])
        if len(cell_counts) < gex_matrices[i][0].shape[1]:
            cell_counts = np.pad(cell_counts, (0, gex_matrices[i][0].shape[1] - len(cell_counts)), 'constant')
            feature_counts = np.pad(feature_counts, (0, gex_matrices[i][0].shape[1] - len(feature_counts)), 'constant')
        tp.append(cell_counts.mean())
        mean_fcs.append(feature_counts.mean())
        med_cs.append(np.median(cell_counts))
        med_fcs.append(np.median(feature_counts))
        gm[gm > 5] = 5
        threshold_mean_cs.append(np.sum(gm, axis=0).mean())

df = pd.DataFrame({
    "tp": tp, "fp": fp, "tech": tech, "ncells": ncells,
    "color": color, "mean_feature_counts": mean_fcs,
    "median_counts": med_cs, "median_feature_counts": med_fcs,
    "threshold_mean_counts": threshold_mean_cs
})
df["total_assigned_counts"] = df["ncells"] * df["tp"]
df.to_csv(args[-6])

# Mean counts plot
plt.figure(figsize=(18, 16))
ax = sns.lineplot(
    x = "fp", y = "tp", hue = "tech", palette = palette, markers = True,
    data = df, marker = "o", markersize = 12, linewidth = 6, legend = False)
plt.ylabel("Avg. counts / cell", fontsize=30)
plt.xlabel("Mutually exclusive coexpression", fontsize=30)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.title("Counts at varying molecule assignment confidence", fontsize=30)
save_path = args[-5]
plt.savefig(save_path, dpi=500)

# Median counts plot
plt.figure(figsize=(18, 16))
ax = sns.lineplot(
    x = "fp", y = "median_counts", hue = "tech", palette = palette, markers = True,
    data = df, marker = "o", markersize = 12, linewidth = 6, legend = False)
plt.ylabel("Median counts / cell", fontsize=30)
plt.xlabel("Mutually exclusive coexpression", fontsize=30)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.title("Counts at varying molecule assignment confidence", fontsize=30)
save_path = args[-4]
plt.savefig(save_path, dpi=500)

# Mean feature counts plot
plt.figure(figsize=(18, 16))
ax = sns.lineplot(
    x = "fp", y = "mean_feature_counts", hue = "tech", palette = palette, markers = True,
    data = df, marker = "o", markersize = 12, linewidth = 6, legend = False)
plt.ylabel("Mean feature counts / cell", fontsize=30)
plt.xlabel("Mutually exclusive coexpression", fontsize=30)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.title("Counts at varying molecule assignment confidence", fontsize=30)
save_path = args[-3]
plt.savefig(save_path, dpi=500)

# Median feature counts plot
plt.figure(figsize=(18, 16))
ax = sns.lineplot(
    x = "fp", y = "median_feature_counts", hue = "tech", palette = palette, markers = True,
    data = df, marker = "o", markersize = 12, linewidth = 6, legend = False)
plt.ylabel("Median feature counts / cell", fontsize=30)
plt.xlabel("Mutually exclusive coexpression", fontsize=30)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.title("Counts at varying molecule assignment confidence", fontsize=30)
save_path = args[-2]
plt.savefig(save_path, dpi=500)

# Threshold at 5 mean counts plot
plt.figure(figsize=(18, 16))
ax = sns.lineplot(
    x = "fp", y = "threshold_mean_counts", hue = "tech", palette = palette, markers = True,
    data = df, marker = "o", markersize = 12, linewidth = 6, legend = False)
plt.ylabel("Mean counts / cell (threshold=5)", fontsize=30)
plt.xlabel("Mutually exclusive coexpression", fontsize=30)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.title("Counts at varying molecule assignment confidence", fontsize=30)
save_path = args[-1]
plt.savefig(save_path, dpi=500)

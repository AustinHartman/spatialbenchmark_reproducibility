# Plot the position of molecules in the cortex and thalamus (Slc17a6, Satb2, and Aqp4)

import sys
import numpy as np
import pandas as pd
from numpy import genfromtxt
import matplotlib.image as mpimg
from matplotlib import pyplot as plt
import random
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import h5py

args = sys.argv[1:]
df = pd.read_csv(args[0])
tx = pd.read_csv(args[1])
ann = pd.read_csv(args[2])
segs = pd.read_csv(args[3])
palette = {
    "Slc17a6": "#d55e00",
    "Satb2": "#cc79a7",
    "Lamp5": "#cc79a7",
    "Neurod6": "#cc79a7",
    "Aqp4": "#009e73",
    "Ntsr2": "#009e73",
}
tx['x'] = tx.x_location
tx['y'] = tx.y_location
tx['gene'] = tx.feature_name
ann_ss = ann[["Unnamed: 0", "predicted.class"]]
ann_ss.columns = ["cell", "annotation"]
celltype_color_mapping = {
    "Neurons": "#d55e005c",
    "Immune": "#6666665c",
    "Oligos": "#6666665c",
    "Astrocytes": "#009e735c",
    "Ependymal": "#6666665c",
    "Vascular": "#6666665c",
    "PeripheralGlia": "#6666665c"
}
ann_ss['color'] = ann_ss.annotation.map(celltype_color_mapping)


xmin=int(3440) - 100
xmax=int(3540) + 100
ymin=int(3350) - 100
ymax=int(3450) + 100
tx_ss = tx[(tx.y_location > ymin) & (tx.y_location < ymax) & (tx.x_location > xmin) & (tx.x_location < xmax)]
tx_ss = tx_ss[tx_ss.feature_name.isin(["Slc17a6", "Aqp4", "Ntsr2"])]
tx_ss["color"] = tx_ss.feature_name.map(palette)
tx_ss.head()
segs_tmp = segs[(segs.vertex_y > ymin-100) & (segs.vertex_y < ymax+100) & (segs.vertex_x > xmin-100) & (segs.vertex_x < xmax+100)]

ann_tmp = ann_ss[ann_ss.cell.isin(segs_tmp.cell_id.unique())]
m = dict(zip(ann_tmp.cell, ann_tmp.color))

patches, colors = [], []
for c in np.unique(segs_tmp.cell_id):
    tmp = segs_tmp[segs_tmp.cell_id == c]
    polygon = Polygon(list(zip(tmp.vertex_x, tmp.vertex_y)),closed=False)
    colors.append(m[c])
    patches.append(polygon)

collection = PatchCollection(patches)
plt.rcParams["figure.figsize"] = (16, 16)
fig,ax = plt.subplots(1)

collection.set_linewidth(3)
collection.set_facecolor(colors)
collection.set_edgecolor([(1, 1, 1, 1) for _ in range(len(patches))])

ax.add_collection(collection)
ax.scatter(tx_ss.x_location, tx_ss.y_location, s=40, alpha=1, c=tx_ss.color)
ax.set_xticks([])
ax.set_yticks([])
ax.set_facecolor("black")
plt.xlim([xmin, xmax])
plt.ylim([ymin, ymax])
plt.savefig(args[-2], dpi=500)

celltype_color_mapping = {
    "Neurons": "#cc79a75c",
    "Immune": "#6666665c",
    "Oligos": "#6666665c",
    "Astrocytes": "#009e735c",
    "Ependymal": "#6666665c",
    "Vascular": "#6666665c",
    "PeripheralGlia": "#6666665c"
}
ann_ss['color'] = ann_ss.annotation.map(celltype_color_mapping)

xmin=int(3470) - 200
xmax=int(3570) - 50
ymin=int(5700) + 400
ymax=int(5800) + 550
tx_ss = tx[(tx.y_location > ymin) & (tx.y_location < ymax) & (tx.x_location > xmin) & (tx.x_location < xmax)]
tx_ss = tx_ss[tx_ss.feature_name.isin(["Satb2", "Neurod6", "Lamp5", "Aqp4", "Ntsr2"])]
tx_ss["color"] = tx_ss.feature_name.map(palette)
tx_ss.head()
segs_tmp = segs[(segs.vertex_y > ymin-100) & (segs.vertex_y < ymax+100) & (segs.vertex_x > xmin-100) & (segs.vertex_x < xmax+100)]

ann_tmp = ann_ss[ann_ss.cell.isin(segs_tmp.cell_id.unique())]
m = dict(zip(ann_tmp.cell, ann_tmp.color))

patches, colors = [], []
for c in np.unique(segs_tmp.cell_id):
    tmp = segs_tmp[segs_tmp.cell_id == c]
    polygon = Polygon(list(zip(tmp.vertex_x, tmp.vertex_y)),closed=False)
    colors.append(m[c])
    patches.append(polygon)

collection = PatchCollection(patches)
plt.rcParams["figure.figsize"] = (16, 16)
fig,ax = plt.subplots(1)

collection.set_linewidth(3)
collection.set_facecolor(colors)
collection.set_edgecolor([(1, 1, 1, 1) for _ in range(len(patches))])

ax.add_collection(collection)
ax.scatter(tx_ss.x_location, tx_ss.y_location, s=40, alpha=1, c=tx_ss.color)
ax.set_xticks([])
ax.set_yticks([])
ax.set_facecolor("black")
plt.xlim([xmin, xmax])
plt.ylim([ymin, ymax])
plt.savefig(args[-1], dpi=500)

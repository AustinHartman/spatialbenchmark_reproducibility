import sys
import numpy as np
import pandas as pd
from numpy import genfromtxt
from matplotlib import pyplot as plt
from PIL import Image, ImageDraw
import random
import csv
import os
import h5py
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
from matplotlib.patches import PathPatch

args = sys.argv[1:]
print(args)

if (args[1][:6] == "vizgen"):
    tx = pd.read_csv(args[0])
    tx = tx[["global_x", "global_y", "gene"]]
    tx.columns = ["x", "y", "gene"] 
elif (args[1][:4] == "tenx"):
    tx = pd.read_csv(args[0])
    tx = tx[["x_location", "y_location", "feature_name", "cell_id"]]
    tx.columns = ["x", "y", "gene", "cell"]
elif (args[1] == "tenx_nuclear"):
    tx = pd.read_csv(args[0])
    tx.loc[tx.overlaps_nucleus == 0, "cell_id"] = -1
    tx = tx[["x_location", "y_location", "feature_name", "cell_id"]]
    tx.columns = ["x", "y", "gene", "cell"]
elif (args[1] == "resolve"):
    # TODO: do molecule assignment
    tx = pd.read_csv(args[0], sep="\t", names=["x", "y", "fov", "gene"])
    tx = tx[["x", "y", "gene"]]
elif (args[1] == "eelfish"):
    tx = pd.read_csv(args[0])
    tx = tx[["c_px_microscope_stitched", "r_px_microscope_stitched", "decoded_genes"]]
    tx.columns = ["x", "y", "gene"]
elif (args[1][:7] == "merfish"):
    tx = pd.read_csv(args[0])
    tx = tx[["global_x", "global_y", "target_gene"]]
    tx.columns = ["x", "y", "gene"]

path = os.path.join("./data/", args[1] + "_molecules.csv")
print(path)
tx.to_csv(path, index=False)

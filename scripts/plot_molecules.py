import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import random

args = sys.argv[1:]
df = pd.read_csv(args[0])
df.loc[df["cell"].isna(), "cell"] = '-1'
df["cell"] = df["cell"].astype(str)
df = df[df['gene'].isin([args[1], args[2]])]
save_path = args[4]
print(save_path)

# plot out all of the molecules in a matplotlib scatter plot with their cell assignments as the color
plt.rcParams["figure.figsize"] = (15, 15)
# create dictionary to map each unique category to a random color
color_dict = {
    args[1]: "#CD2029",
    args[2]: "#3851A1",
}

# map each category to its corresponding color
df['color'] = df['gene'].map(color_dict)
plt.scatter(df.x, df.y, c=df.color, s=float(args[3]))
plt.title("Molecules colored by gene", fontsize=35)
plt.axis('off')
plt.savefig(save_path, dpi=300)

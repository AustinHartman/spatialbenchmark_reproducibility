import sys
import os
import numpy as np
import pandas as pd

args = sys.argv[1:]
print(args)

tx = pd.read_csv(args[0])

if (args[1] == "vizgen" and args[2] == "cortex"):
    tx = tx[(tx.x > 3500) & (tx.x < 4500) & (tx.y > 500) & (tx.y < 1500)]
elif (args[1] == "vizgen_rep2" and args[2] == "cortex"):
    tx = tx[(tx.x > 3000) & (tx.x < 4000) & (tx.y > 8250) & (tx.y < 9250)]
elif (args[1] == "vizgen_rep3" and args[2] == "cortex"):
    tx = tx[(tx.x > 5500) & (tx.x < 6500) & (tx.y > 6000) & (tx.y < 7000)]
elif ((args[1] == "tenx" or args[1] == "tenx_nuclear") and args[2] == "cortex"):
    tx = tx[(tx.x > 3300) & (tx.x < 4300) & (tx.y > 5600) & (tx.y < 6600)]
elif (args[1] == "tenx_rep2" and args[2] == "cortex"):
    tx = tx[(tx.x > 7000) & (tx.x < 8000) & (tx.y > 5400) & (tx.y < 6400)]
elif (args[1] == "tenx_rep3" and args[2] == "cortex"):
    tx = tx[(tx.x > 5500) & (tx.x < 6500) & (tx.y > 5850) & (tx.y < 6850)]
elif (args[1] == "resolve" and args[2] == "cortex"):
    tx = tx[(tx.x > 15000) & (tx.x < 20000) & (tx.y > 2000) & (tx.y < 7000)]
elif (args[1] == "eelfish" and args[2] == "cortex"):
    tx = tx[(tx.x > 33000) & (tx.x < 38000) & (tx.y > 30000) & (tx.y < 35000)]
elif (args[1] == "merfish" and args[2] == "cortex"):
    tx = tx[(tx.x > 1400) & (tx.x < 2250) & (tx.y > 4400) & (tx.y < 5250)]
elif (args[1] == "merfish_rep2" and args[2] == "cortex"):
    tx = tx[(tx.x > -3200) & (tx.x < -2000) & (tx.y > 1200) & (tx.y < 2000)]
elif (args[1] == "vizgen" and args[2] == "thalamus"):
    tx = tx[(tx.x > 3250) & (tx.x < 4000) & (tx.y > 3250) & (tx.y < 4000)]
elif ((args[1] == "tenx" or args[1] == "tenx_nuclear") and args[2] == "thalamus"):
    tx = tx[(tx.x > 3100) & (tx.x < 4400) & (tx.y > 2800) & (tx.y < 4100)]
elif (args[1] == "resolve" and args[2] == "thalamus"):
    tx = tx[(tx.x > 20000) & (tx.x < 25000) & (tx.y > 17000) & (tx.y < 22000)]
elif (args[1] == "eelfish" and args[2] == "thalamus"):
    tx = tx[(tx.x > 42000) & (tx.x < 47000) & (tx.y > 18000) & (tx.y < 23000)]
elif (args[1] == "merfish" and args[2] == "thalamus"):
    tx = tx[(tx.x > -1000) & (tx.x < 0) & (tx.y > 5500) & (tx.y < 6500)]
elif (args[1] == "de_tenx" and args[2] == "cortex"):
    tx = tx[(tx.x > 3000) & (tx.x < 7500) & (tx.y > 5700) & (tx.y < 6500)]
elif (args[1] == "de_tenx" and args[2] == "thalamus"):
    tx = tx[(tx.x > 3100) & (tx.x < 6700) & (tx.y > 2400) & (tx.y < 3900)]

path = os.path.join("./data/", args[1] + "_" + args[2] + "_molecules.csv")
tx.to_csv(path, index=False)

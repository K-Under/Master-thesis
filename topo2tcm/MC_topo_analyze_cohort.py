import MC_topo_analyzer as mca
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# importing auxillary scripts
import harry_plotter as hp
import geometry_funcs as gf


"""
Cohort analysis of recon files using MC_topo_analyzer.py
"""


path = Path(r"D:\projection_data_sample\Recon\2024-10-01_13-45_L_Standard_segm_0.2_Result")

data = mca.load_sim(path)

# Access image data
print(len(data))

hist1, hist2 = gf.slyce2histogram(data[750])

hp.showim(data)
hp.plot_data_simple(hist1)
hp.plot_data_simple(hist2)

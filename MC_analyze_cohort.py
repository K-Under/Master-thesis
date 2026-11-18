import MC_analyzer as mca
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

"""
Cohort analysis of .mifx and .raw files using MC_analyzer.py
"""
#NOTETOSELF: 0 = AIR, 9 = RING, CHECK .bat FILE FOR MORE INFO


data_path = Path("D:\InputData")
for f in data_path.iterdir():
    if not ("Standard_segm" in f.name and f.suffix) == ".mifx":
        continue
    print(f"Loading {f}...")
    data = mca.load_mifx(f)
    breast = mca.BreastAnalyzer(data)
    # breast.plot_data_simple(breast.diameter_lyst)
    # breast.plot_data_simple(breast.dose_volume_reduction_factor())
    #   tcm_path = f.with_suffix(".tcm")
    # breast.tcm_file_write(tcm_path, "testdir")
    #   breast.tcm_file_geometry(tcm_path, "step_tcm", geom_bool=False)
    breast.showim(data)
print("Done!")





 # testing the breast class

# mifx_path = Path(r"D:\InputData\2024-10-01_13-45_L_Standard_segm.mifx")
# data = mca.load_mifx(mifx_path)
# breast = mca.BreastAnalyzer(data)
#
# tcm_path = Path(r"C:\Users\Pepponita\Documents\Uni\Masterarbeit\Test\test2.tcm")
# breast.tcm_file_write(tcm_path)
#
# red_fac = breast.dose_volume_reduction_factor()
# dose_est = breast.dose_volume_estimate()
#
# x_axis = np.linspace(0,len(red_fac),len(red_fac))
#
# plt.plot(x_axis,red_fac)
# plt.plot(x_axis,dose_est)
# plt.show()






# old mean glandular dose shower

# breast_dose = mca.ring_remover(data)
# MGD = mca.dose_vol_estimate(breast_dose)
# x_axis = np.linspace(0,len(MGD),len(MGD))
#
# plt.plot(x_axis,MGD)
# plt.show()



# used to check the shin paper (MGD multiple sims)

# length_lyst = []
# mean_dose_lyst = []
# data_path = Path("D:\Sim04_Segmented_Modulated_Loop\MC")
# for f in data_path.iterdir():
#     if "Standard_segm" in f.name and f.suffix == ".mifx":
#         data = mca.load_mifx(f)
#         breast = mca.BreastAnalyzer(data)
#         length_lyst.append(breast.breast_length)
#         breast_dose = breast.breast_data
#         mean_dose_lyst.append(breast.dose_mean(breast_dose, True))
#
# plt.scatter(length_lyst, mean_dose_lyst)
# plt.xlabel("Breast length [mm]")
# plt.ylabel('Mean absorbed dose (MAD) [mGy]')
# plt.title("MAD & Breast length correlation")
# plt.grid()
# plt.show()



# used for length and width estimation

# mifx_path = Path(r"D:\InputData\2024-10-01_13-45_L_Standard_segm.mifx")
# data = mca.load_mifx(mifx_path)
#
# length = mca.breast_length_estimate(data)
# print(length)
#
# mca.show_width_plot(data)

# width = mca.breast_slice_width_estimate(data[300])
# print(f"Breast diameter: {width} mm or {width*10/3} pixels")

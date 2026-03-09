import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns

# Read averaged csv file
dfLA = pd.read_csv(r"C:/Users/hwest10/ICU-simulator/manuscript/Mean_abandonment_array_LA_853333.csv")
dfLB = pd.read_csv(r"C:/Users/hwest10/ICU-simulator/manuscript/Mean_abandonment_array_LB_342595.csv")
dfLC = pd.read_csv(r"C:/Users/hwest10/ICU-simulator/manuscript/Mean_abandonment_array_LC237246.csv")

# Get dimensions of array
n_hcp_L = len(dfLA.columns)
n_beds_L = dfLA.shape[0]

# Get appropriately indexed lists
xtix_L = np.linspace(1, n_hcp_L, n_hcp_L, dtype=int)
ytix_L = np.linspace(1, n_beds_L, n_beds_L, dtype=int)

# Get min and max values to normalize heatmaps
vmin = min(dfLA.min(axis=None),dfLB.min(axis=None),dfLC.min(axis=None))
vmax = max(dfLA.max(axis=None),dfLB.max(axis=None),dfLC.max(axis=None))

# Plot codes----------------------------------------------------------------
# open subplot figure with 3x1 size High resource or Low resource
fig, ax = plt.subplots(nrows=1,ncols=3, sharex=False, layout='constrained')
fig.suptitle('Low resource facility: Accumulated queue abandonment over 365 days')
fig.supxlabel('Clinicians on staff')
fig.supylabel('ICU Beds')

# Set parameters for plot
plt.rcParams["figure.figsize"] = (18,8)

# Set universal color range
#cmap = mcolors.LinearSegmentedColormap.from_list("n",["#B3BCFC","#5E6DDF","#4654C4","#2C3AA8","#1B288B","#121F7C",
#                       "#081260","#040D55","#00094A","#000000"])

cmap = mcolors.LinearSegmentedColormap.from_list("n",["#00876c","#45a074","#72b97c","#9fd184","#cee98f","#ffff9d",
                       "#fedb79","#fab560","#f38f52","#e7674e","#d43d51"])
norm = plt.Normalize(-100,100)

# Create heatmap
s1 = sns.heatmap(dfLA, cmap=cmap,fmt=".1f", vmin=vmin, vmax=vmax, ax=ax[0], annot=False, xticklabels=xtix_L, yticklabels=ytix_L, cbar=False)
s2 = sns.heatmap(dfLB, cmap=cmap,fmt=".1f", vmin=vmin, vmax=vmax, ax=ax[1], annot=False, xticklabels=xtix_L, yticklabels=ytix_L, cbar=False)
s3 = sns.heatmap(dfLC, cmap=cmap,fmt=".1f", vmin=vmin, vmax=vmax, ax=ax[2], annot=False, xticklabels=xtix_L, yticklabels=ytix_L)

# Set titles
ax[0].set_title('Scenario A: Pre-epidemic')
ax[1].set_title('Scenario B: Initial outbreak')
ax[2].set_title('Scenario C: Increased transmissibility')

# Invert y axes
ax[0].invert_yaxis()
ax[1].invert_yaxis()
ax[2].invert_yaxis()

# Make axis labels smaller
ax[0].tick_params(axis='both', which='major', labelsize=7)
ax[1].tick_params(axis='both', which='major', labelsize=7)
ax[2].tick_params(axis='both', which='major', labelsize=7)

lin_x = [1,2,3,4,5,6,7,8]
lin_y = [18 - 2*x for x in lin_x]

sns.lineplot(x=lin_x,y=lin_y,ax=ax[0],color='black',drawstyle='steps-pre',label='Budget threshold')
sns.lineplot(x=lin_x,y=lin_y,ax=ax[1],color='black',drawstyle='steps-pre',label='Budget threshold')
sns.lineplot(x=lin_x,y=lin_y,ax=ax[2],color='black',drawstyle='steps-pre',label='Budget threshold')

plt.legend()

# Show the plot
plt.show()
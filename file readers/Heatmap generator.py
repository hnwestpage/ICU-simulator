import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read csv file
df = pd.read_csv(r"C:\Users\hwest10\ICU-simulator\Mean_abandonment_array.csv")

# Set parameters for plot
plt.rcParams["figure.figsize"] = (18,8)

title = 'Test output budget array'
x_axis = 'Clinicians'
y_axis = 'Beds'

plt.title(title,fontsize=14)

# Create heatmap
sns.heatmap(df, cmap="rocket_r",fmt=".1f")
plt.show()
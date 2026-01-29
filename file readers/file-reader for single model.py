# File reader to read results off of cluster
import numpy as np
import math
import matplotlib.pyplot as plt

# Lists
t_ints = list(range(0,365,1))       # Averaging function time intervals
run_t = []
run_X = []
run_C = []
run_Z = []
run_F = []
run_H = []
run_J = []
run_q_D = []
run_q_D_F = []
run_icu_D = []
run_N_T = []
R = 0

# Functions
def find_avg_run(t,X):
  avg_X = []

  for time in t_ints:

    X_vals = []

    for row in t:

      k = 0
      r_point = 0

      while row[k] <= time:
        r_point = row[k]
        k += 1

      x_index = t.index(row)
      y_index = row.index(r_point)
      X_vals.append(X[x_index][y_index])

    avg = np.mean(X_vals)
    avg_X.append(avg)

  return avg_X

# Create multipanel figure
fig, ax = plt.subplots(nrows=4,ncols=1, sharex=True, layout='constrained')

# Opening print file
f = open('slurm-5198311.out')

# Reading each line of printed file
for line in f:
    if 'Run number' in line:
        run_str = line.split(' ')[2]
        run = int(run_str)
        print(run)

    if 'Timestamp list:' in line:
        line_clean = line.replace('[', '').replace(']','')
        t_str = line_clean.split(' ')[2:]
        t = [float(str(i).replace(",", "").replace("\n","")) for i in t_str]
        run_t.append(t)

    if 'Number in queue' in line:
        line_clean = line.replace('[', '').replace(']','')
        X_str = line_clean.split(' ')[3:]
        X = [float(str(i).replace(",", "").replace("\n","")) for i in X_str]
        run_X.append(X)
        ax[0].plot(t,X, color='lightgray')

    if 'Number of infectious patients in queue' in line:
        line_clean = line.replace('[', '').replace(']','')
        C_str = line_clean.split(' ')[6:]
        C = [float(str(i).replace(",", "").replace("\n","")) for i in C_str]
        run_C.append(C)
        ax[0].plot(t,C, color='lightpink')

    if 'Number in ICU' in line:
        line_clean = line.replace('[', '').replace(']','')
        Z_str = line_clean.split(' ')[3:]
        Z = [float(str(i).replace(",", "").replace("\n","")) for i in Z_str]
        run_Z.append(Z)
        ax[1].plot(t,Z, color='lightgray')

    if 'Number of infectious patients in ICU' in line:
        line_clean = line.replace('[', '').replace(']','')
        F_str = line_clean.split(' ')[6:]
        F = [float(str(i).replace(",", "").replace("\n","")) for i in F_str]
        run_F.append(F)
        ax[1].plot(t,F, color='palegreen')

    if 'Number of active clinicians' in line:
        line_clean = line.replace('[', '').replace(']','')
        H_str = line_clean.split(' ')[4:]
        H = [float(str(i).replace(",", "").replace("\n","")) for i in H_str]
        run_H.append(H)
        ax[3].plot(t,H, color='lightgray')

    if 'Queue departures' in line:
        line_clean = line.replace('[', '').replace(']','')
        q_D_str = line_clean.split(' ')[2:]
        q_D = [float(str(i).replace(",", "").replace("\n","")) for i in q_D_str]
        run_q_D.append(q_D)
        ax[2].plot(t,q_D, color='lightgray')

    if 'Infectious queue departures' in line:
        line_clean = line.replace('[', '').replace(']','')
        q_D_F_str = line_clean.split(' ')[3:]
        q_D_F = [float(str(i).replace(",", "").replace("\n","")) for i in q_D_F_str]
        run_q_D_F.append(q_D_F)
        ax[2].plot(t,q_D_F, color='lightblue')



# Calculating the average run
X_avg = find_avg_run(run_t, run_X)
C_avg = find_avg_run(run_t, run_C)
Z_avg = find_avg_run(run_t, run_Z)
F_avg = find_avg_run(run_t, run_F)
H_avg = find_avg_run(run_t, run_H)
q_D_avg = find_avg_run(run_t, run_q_D)
q_D_F_avg = find_avg_run(run_t, run_q_D_F)
#N_T_avg = find_avg_run(run_t, run_N_T)

ax[0].plot(t_ints, X_avg, color='red', linewidth = 2)
ax[0].plot(t_ints, C_avg, color='lightcoral', linewidth = 2)
ax[1].plot(t_ints, Z_avg, color='green', linewidth = 2)
ax[1].plot(t_ints, F_avg, color='limegreen', linewidth = 2)
#ax[1].plot(t_ints, N_T_avg, 'k:', linewidth = 2)
ax[2].plot(t_ints, q_D_avg, color='black', linewidth = 2)
ax[2].plot(t_ints, q_D_F_avg, color='dimgray', linewidth = 2)
ax[3].plot(t_ints, H_avg, color='blue', linewidth = 2)

fig.suptitle('M/M/'+str(30)+' with dynamic arrival rate, HCP susceptible')
fig.supxlabel('Time')

ax[0].set_ylabel('Patients in queue')
ax[1].set_ylabel('Patients in ICU')
ax[2].set_ylabel('Queue Departures')
ax[3].set_ylabel('Active Clinicians')

#fig.tight_layout()
fig.set_size_inches(7,9)
plt.show()
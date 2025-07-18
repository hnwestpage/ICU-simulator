# ICU PREPAREDNESS SIMULATOR 2.1-------------------------------------

# imports
import numpy as np
import math
#---------------------------------------------------------------------
# PRELIMINARY FUNCTIONS
# Budget operations calculator

# Function to find linear combinations of HCP and beds given a fixed budget
def get_budget_options(U, a, b):
  # Get max values for budget
  HCP_max = int(U/a)
  n_max = int(U/b)

  budget_array = []
  # Build array of pairs
  for i in range(n_max - 1):
    inner_vec = []
    for j in range(HCP_max - 1):
      pair = (i+1,j+1)
      inner_vec.append(pair)
    budget_array.append(inner_vec)

  return budget_array

# Parameter initializations

# Resources                # Description
U = 2.6 * 10**7            # Hospital budget in USD (total expense = $26 million)
a = 3.0 * 10**5            # Average HCP yearly salary ($300,000)
b = 1.25 * 10**5           # Average cost of one ICU bed for a year ($125,000)
c = 0.35                   # Percent of budget spent on constant costs (Utilities, maintenence)
d = 0.27                   # Percent of budget spent on nursing personnel

# Get combinations based on allotted budget
budget_array = get_budget_options((1-c-d)*U, a, b)

# Rates
lmbda = 7            # Baseline arrival rate of patients to ICU
mu_1 = 1/3.3         # Departure rate of baseline (non infectious) patients from ICU (1/recovery)
mu_2 = 1/8.0         # Departure rate of COVID-19 (infectious) patients from ICU (1/recovery)
rho_1 = 1/20         # Departure rate of baseline patients from the queue (renege)
rho_2 = 1/20         # Departure rate of COVID-19 patients from the queue (renege)
eta = 0.020          # Rate at which HCP becomes infected in the workplace
nu = 1/6             # Rate at which infected HCP recover/return to work

# Proportions
mort_ICU = 0.245           # Proportion of baseline patients who perish in the ICU- Kim et al. (2021), Vincent et al. (2020)
mort_Q = 0.507             # Proportion of baseline patients who perish in the queue- Boumendil et al. (2012)
mort_ICU_19 = 0.245        # Proportion of COVID-19 patients who perish in the ICU- !NEEDS REF
mort_Q_19 = 0.507          # Proportion of COVID-19 patients who perish in the queue- !NEEDS REF
crit = 0.00102             # Proportion of infected indiviuals requiring ICU- Moon et al., (2022)
r =  1/9.5                 # Threshold ratio of HCP to beds (servers)

# SEIR parameters - Original Wuhan
beta = 0.0000033    # contact rate between susceptibles and infectives
gamma = 1/5.2       # reciprocal of mean exposed period (5.2 days)
alpha = 1/8         # reciprocal of mean recovery period (8 days)
tt = 6              # number of time steps per day (~4 hour time steps)

# Miscellaneous
N = 1*10**5                       # population size (urban setting baseline)
M = 1                            # Number of sample paths
T = 365                           # Time (days)
t_ints = list(range(0,T,1))       # Averaging function time intervals
x = np.linspace(0, T, T*tt + 1)   # SEIR timescale variable

# Set and Get Functions

# Function to check the current number of infected individuals
def check_I(current_t):
  l_point = 0

  for time in x:
    if time <= current_t:
      l_point = time
    else:
      break

  index = np.where(x == time)[0][0]
  num = I[index]

  return num

# Function to check the current arrival rate
def check_lambda(t, lmbda):
  return lmbda + crit*check_I(t)

# Function to check the current ICU departure rates
#   i = total # patients in ICU, j = # infectious
def check_mu(i, j, n_current, mu_1, mu_2):
  if i <= 0:
    return 0
  else:             # revisit this choice (what if we are over capacity?)
    mu_1_t = (i-j)*mu_1
    mu_2_t = j*mu_2

  return mu_1_t + mu_2_t


# Function to check the current queue departure rate
#   i = total patients in queue (X_occ length)
#   j = infectious patients in queue (X_occ sum)
def check_rho(X_occ, rho_1, rho_2):
  i = len(X_occ)
  j = sum(X_occ)
  if i+j <= 0:
    return 0
  else:
    rho_1_t = i*rho_1
    rho_2_t = j*rho_2

  return rho_1_t + rho_2_t

# Function to check the current HCP infection rate
#   i = number of infectious patients
#   H = number of clinicians on staff
def check_eta(i, H, eta):
  if i > 0 and H > 0:
    return i*eta
  else:
    return 0

# Function to check the current HCP return rate
#   i = number of infected HCP
def check_nu(i, nu):
  if i <= 0:
    return 0
  else:
    return i*nu

# Function to check current unit capacity based on active HCP
def check_capacity(n, H, N_T, r):
  N_new = np.floor(H[-1]*(1/r))
  N_update = np.min((N_new, n))
  N_T.append(N_update)

  return N_T

# Event Functions

# Simulated SEIR-pandemic
def SEIR(S_0, E_0, I_0, R_0, days, N): #Inputs are initial no. of patients and days to run

  #Calculating basic reproduction number
  R0 = round((N*beta)/(alpha),2)

  # Initialize lists to record number of patients in compartment at each time step
  S = [S_0]
  E = [E_0]
  I = [I_0]
  R = [R_0]

  # Begin for loop to iterate through each time step for duration of model
  for t in range(days*tt):

    # Implementing previous equations
    dS = (-beta*S[t]*I[t] + (alpha)*I[t]) / tt
    dE = (beta*S[t]*I[t] - gamma*E[t]) / tt
    dI = (gamma*E[t] -(alpha)*I[t]) / tt
    dR = (alpha)*I[t] / tt

    # Adding the condition that values must be non-negative to be recorded
    if E[t] + dE >= 0:
      E.append(E[t] + dE)
    else:
      E.append(0)

    if I[t] + dI >= 0:
      I.append(I[t] + dI)
    else:
      I.append(0)

    R.append(R[t] + dR)

    # N = S + E + I + R
    S.append(N - E[t] - I[t] - R[t])

  # Save x as list of uniform timestamps
  x = np.linspace(0, days, days*tt + 1)

  return I

# Event 1: An arrival to the queue ***UPDATE COMPLETE
def arrival(t, n_current, X_occ, Z, F, H, J, q_D, q_D_F, icu_D):
  # Determine patient type------------------------------------------------------
  # Summing rate for infectious and non-infectious arrival
  non_inf = lmbda
  inf = crit*check_I(t)
  rate = non_inf + inf

  # Selection of patient type: 0 = non-infectious, 1 = infectious
  events = [0,1]
  probs = [non_inf/rate, inf/rate]
  y = np.random.choice(a=events, p=probs)
  #y = y.astype(int)

  # Determine where to send arrival
  if Z[-1] >= n_current:              # if all beds are at capacity
    X_occ.append(int(y))                   # the new arrival waits in the queue
    Z.append(int(Z[-1]))
    F.append(int(F[-1]))
    H.append(int(H[-1]))
    J.append(int(J[-1]))
    q_D.append(int(q_D[-1]))
    q_D_F.append(int(q_D_F[-1]))
    icu_D.append(int(icu_D[-1]))
  else:                               # otherwise there is an empty bed
    #X.append(int(X[-1]))
    Z.append(int(Z[-1]+1))            # add the new arrival to the open bed
    F.append(int(F[-1]+y))            # add value of rv to update infectious counts
    H.append(int(H[-1]))
    J.append(int(J[-1]))
    q_D.append(int(q_D[-1]))
    q_D_F.append(int(q_D_F[-1]))
    icu_D.append(int(icu_D[-1]))
  return X_occ, Z, F, H, J, q_D, q_D_F, icu_D

# Event 2: A departure from the queue (renege)  ***UPDATE COMPLETE
def q_depart(X_occ, Z, F, H, J, q_D, q_D_F, icu_D):
  # Determine patient type
  non_inf = len(X_occ) - sum(X_occ)
  inf = sum(X_occ)
  rate = non_inf + inf

  # Selection of patient type: 0 = non-infectious, 1 = infectious
  events = [0,1]
  probs = [non_inf/rate, inf/rate]
  y = np.random.choice(a=events, p=probs)

  if X_occ:                   # if there was someoine in the queue to begin with
    #X.append(int(X[-1]-1))        # reduce number in queue by 1
    X_occ.remove(y)                # removes first occurence of rv
    Z.append(int(Z[-1]))
    F.append(int(F[-1]))
    H.append(int(H[-1]))
    J.append(int(J[-1]))
    q_D.append(int(q_D[-1]+1))      # record one queue departure
    q_D_F.append(int(q_D_F[-1]+y))  # record patient type
    icu_D.append(int(icu_D[-1]))
  else:
    print("Something went wrong: the queue has negative people.")

  return X_occ, Z, F, H, J, q_D, q_D_F, icu_D

# Event 3: A departure from the ICU***************UPDATE COMPLETE
def icu_depart(n_current, X_occ, Z, F, H, J, q_D, q_D_F, icu_D):
  # Determine departing patient type
  non_inf = Z[-1] - F[-1]
  inf = F[-1]
  rate = non_inf + inf

  # Selection of departing patient type: 0 = non-infectious, -1 = infectious
  events = [0,-1]
  probs = [non_inf/rate, inf/rate]
  y = np.random.choice(a=events, p=probs)

  # potentially appending something here, but maybe not

  if X_occ and Z[-1] <= n_current:  # if there was someone waiting in the queue & max is not exceeded
    z = X_occ[0]               # next patient(value either 0 or 1)
    X_occ.remove(z)            # move the next person from the queue to the empty bed
    Z.append(int(Z[-1]))       # total number in ICU remains the same
    F.append(int(F[-1]+y+z))     # number of infectious patients updated based on next patient
    H.append(int(H[-1]))
    J.append(int(J[-1]))
    q_D.append(int(q_D[-1]))
    q_D_F.append(int(q_D_F[-1]))
    icu_D.append(int(icu_D[-1]+1))    # record one ICU departure
  else:                               # otherwise, if no one is in the queue
    Z.append(int(Z[-1]-1))            # reduce the number in the ICU by 1
    F.append(int(F[-1]+y))            # update no. infectious by depating patient type
    H.append(int(H[-1]))
    J.append(int(J[-1]))
    q_D.append(int(q_D[-1]))
    q_D_F.append(int(q_D_F[-1]))
    icu_D.append(int(icu_D[-1]+1))    # record one ICU departure
  return X_occ, Z, F, H, J, q_D, q_D_F, icu_D

# Event 4: A HCP becomes infected
def HCP_infect(X_occ, Z, F, H, J, q_D, q_D_F, icu_D):
  #X.append(int(X[-1]))
  Z.append(int(Z[-1]))
  F.append(int(F[-1]))
  H.append(int(H[-1]-1))        # reduce the number of active HCP by 1
  J.append(int(J[-1]+1))        # increase the number of infected HCP by 1
  q_D.append(int(q_D[-1]))
  q_D_F.append(int(q_D_F[-1]))
  icu_D.append(int(icu_D[-1]))
  return X_occ, Z, F, H, J, q_D, q_D_F, icu_D

# Event 5: A HCP recovers and returns to work
def HCP_return(X_occ, Z, F, H, J, q_D, q_D_F, icu_D):
  #X.append(int(X[-1]))
  Z.append(int(Z[-1]))
  F.append(int(F[-1]))
  H.append(int(H[-1]+1))        # increase the number of active HCP by 1
  J.append(int(J[-1]-1))        # reduce the number of infected HCP by 1
  q_D.append(int(q_D[-1]))
  q_D_F.append(int(q_D_F[-1]))
  icu_D.append(int(icu_D[-1]))
  return X_occ, Z, F, H, J, q_D, q_D_F, icu_D

# Finding times and averages

# Function which finds the average value of each randomized run at uniform timesteps in [0,T].
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

# Function to find the peak of the infection wave
def find_peak(x, I):
  max = -1
  for i in I:
    if i > max:
      max = i
    else:
      max = max

  max_index = I.index(max)
  t_max = x[max_index]
  return t_max

# function to average the number of queue departures/deaths
def get_mortality_stats(run_q_D):
  entries = []                          # empty list for grabbing last entry of run

  for i in range(M):                    # for each run
    entries.append(run_q_D[i][-1])      # append the last entry of the run

  avg = np.mean(entries)
  #print('Queue departure average: '+str(avg))

  std = np.std(entries)
  #print('Proportion of mortality: '+str(mort_Q))

  U_mort_mean = mort_Q * avg
  #print('Reported untreated deaths: '+str(U_mort_mean))

  U_mort_std = mort_Q * std

  return U_mort_mean, U_mort_std

# Function to find the (beds/hcp combo) minimizer of queue deaths
def find_mort_min(U_mort_avg):
  x_star = np.argmin(U_mort_avg)
  return x_star

#---------------------------------------------------------------------
# SIMULATOR LOOP

# SEIR model
# Seeding S, E, I, R
E_0 = 0.001*N      # 0.1% of population exposed
S_0 = N - E_0
I_0 = 0
R_0 = 0

# Implmenting model: S_0, E_0, I_0, days
I = SEIR(S_0, E_0, I_0, R_0, T, N)

# SINGLE SCENARIO MODEL****************************************************************************************************************
def single_model(H, n, lmbda, mu_1, mu_2, rho_1, rho_2, eta, nu, ratio):

  # Single scenario cell- fixed HCP and n
  H_current = H
  n_current = n

  # Print baseline traffic density lambda/(n x mu)
  traffic = lmbda/(n_current*mu_1)
  print('The baseline traffic density is '+str(traffic)+'%')

  # Initializing empty arrays for tracking all sample runs
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

  for i in range (M):
    # Seeding initial conditions
    # Queue placeholder
    X_occ = []

    # Initial condition values
    t_0 = 0                                                   # Current time
    X_0 = len(X_occ)                                          # Number in queue
    C_0 = sum(X_occ)                                          # Number COVID-19 patients in queue
    Z_0 = np.min((math.ceil(traffic*n_current), n_current))   # Number in ICU (seeding at expected traffic density occupancy)
    F_0 = 0                                                   # Number of infectious (COVID-19) patients in ICU
    H_0 = H_current                                           # Number of active HCP
    J_0 = 0                                                   # Nubmber of inactive HCP
    q_D_0 = 0                                                 # Number of queue departures (untreated)
    q_D_F_0 = 0                                               # Number of queue departures who were infectious
    icu_D_0 = 0                                               # Number of ICU departures (treated)
    n_0 = n_current                                           # Current ICU size maximum

    # Queue placeholder
    X_occ = []

    # Initializing list forms with initial conditions
    t = [t_0]                 # time stamp list
    X = [X_0]                 # Number in queue at time t
    C = [C_0]                 # Number COVID-19 patients in queue at time t
    Z = [Z_0]                 # Number of total patients in ICU at time t
    F = [F_0]                 # Number of infectious patients i ICU at time t
    H = [H_0]                 # Number of active HCP at time t
    J = [J_0]                 # Number of inactive HCP at time t
    q_D = [q_D_0]             # Number of queue departures at time t (untreated)
    q_D_F = [q_D_F_0]         # Number of queue departures who were infectious
    icu_D = [icu_D_0]         # Number of ICU departures at time t (treated)
    N_T = [n_0]               # Current ICU size maximum at time t


    # While the current time is less than the end time
    while t[-1] < T:

      # Check capacity of the unit based on active HCP
      N_T = check_capacity(n_current, H, N_T, ratio)

      # Find the time that the next event occurs
      # Check all rates
      lmbda_t = check_lambda(t[-1], lmbda)
      mu_t = check_mu(Z[-1],F[-1],N_T[-1], mu_1, mu_2)
      rho_t = check_rho(X_occ, rho_1, rho_2)
      eta_t = check_eta(F[-1],H[-1], eta)
      nu_t = check_nu(J[-1], nu)

      # Total rate
      rate = lmbda_t + mu_t + rho_t + eta_t + nu_t

      # Time update
      delta = np.random.exponential(scale=1/rate,size=1)
      t.append(float(t[-1] + delta))

      # Determine which event occurs at updated time
      events = [1, 2, 3, 4, 5]
      probs = [lmbda_t/rate, rho_t/rate, mu_t/rate, eta_t/rate, nu_t/rate]
      Y = np.random.choice(a=events, p=probs)

      # Event occurence
      if Y == 1:
        X_occ, Z, F, H, J, q_D, q_D_F, icu_D = arrival(t[-1], N_T[-1], X_occ, Z, F, H, J, q_D, q_D_F, icu_D)
      elif Y == 2:
        X_occ, Z, F, H, J, q_D, q_D_F, icu_D = q_depart(X_occ, Z, F, H, J, q_D, q_D_F, icu_D)
      elif Y == 3:
        X_occ, Z, F, H, J, q_D, q_D_F, icu_D = icu_depart(N_T[-1], X_occ, Z, F, H, J, q_D, q_D_F, icu_D)
      elif Y == 4:
        X_occ, Z, F, H, J, q_D, q_D_F, icu_D = HCP_infect(X_occ, Z, F, H, J, q_D, q_D_F, icu_D)
      elif Y == 5:
        X_occ, Z, F, H, J, q_D, q_D_F, icu_D = HCP_return(X_occ, Z, F, H, J, q_D, q_D_F, icu_D)
      else:
        print("Something went wrong. :c")

      # Determine length of queue
      X.append(len(X_occ))

      # Determine number of infectious patients in queue
      C.append(sum(X_occ))

    # Record run
    run_t.append(t)
    run_X.append(X)
    run_C.append(C)
    run_Z.append(Z)
    run_F.append(F)
    run_H.append(H)
    run_q_D.append(q_D)
    run_q_D_F.append(q_D_F)
    run_N_T.append(N_T)

  # Printing resulting array of runs
  print('Timestamp list:',run_t)
  print('Number in queue',X)
  print('Number of infectious patients in queue',C)
  print('Number in ICU',Z)
  print('Number of infectious patients in ICU',F)
  print('Number of active clinicians',H)
  print('Queue departures',q_D)
  print('Infectious queue departures',q_D_F)
  print('Effective ICU Capacity',N_T)

  """ Calculating the average run: commented out for now, may be able to do this locally
  X_avg = find_avg_run(run_t, run_X)
  C_avg = find_avg_run(run_t, run_C)
  Z_avg = find_avg_run(run_t, run_Z)
  F_avg = find_avg_run(run_t, run_F)
  H_avg = find_avg_run(run_t, run_H)
  q_D_avg = find_avg_run(run_t, run_q_D)
  q_D_F_avg = find_avg_run(run_t, run_q_D_F)
  N_T_avg = find_avg_run(run_t, run_N_T)
  """

  # Untreated mortality mean and standard deviation
  mean, stdev = get_mortality_stats(run_q_D)
  U_mort_avg = mean
  U_mort_std = stdev

  return float(U_mort_avg), float(U_mort_std)

# MULTI SCENARIO MODEL*******************************************************************************************************************
def multi_model(U, a, b, mu_1, mu_2, eta, nu):
  # Get combinations based on allotted budget
  budget_array = get_budget_options(U, a, b)

  # Initializing logs of untreated deaths per budget option
  U_mort_avg = np.zeros((int(U/b),int(U/a)))      # dimension (n x H)
  U_mort_std = np.zeros((int(U/b),int(U/a)))    # dimension (n x H)

  # Per budget scenario
  for row in budget_array:
    for pair in row:
      # Establishing the current scenario
      n_current = pair[0]
      H_current = pair[1]

      # Pre-pandemic traffic density
      traffic = lmbda/(n_current*mu_1)

      # Initializing empty arrays for tracking all sample runs
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

      # For each sample run
      for i in range (M):
        # Seeding initial conditions
        # Queue placeholder
        X_occ = []

        # Initial condition values
        t_0 = 0                                                   # Current time
        X_0 = len(X_occ)                                          # Number in queue
        C_0 = sum(X_occ)                                          # Number COVID-19 patients in queue
        Z_0 = np.min((math.ceil(traffic*n_current), n_current))   # Number in ICU (seeding at expected traffic density occupancy)
        F_0 = 0                                                   # Number of infectious (COVID-19) patients in ICU
        H_0 = H_current                                           # Number of active HCP
        J_0 = 0                                                   # Number of inactive HCP
        q_D_0 = 0                                                 # Number of queue departures (untreated)
        q_D_F_0 = 0                                               # Number of queue departures who were infectious
        icu_D_0 = 0                                               # Number of ICU departures (treated)
        n_0 = n_current                                           # Current ICU size maximum

        # Queue placeholder
        X_occ = []

        # Initializing list forms with initial conditions
        t = [t_0]                 # time stamp list
        X = [X_0]                 # Number in queue at time t
        C = [C_0]                 # Number COVID-19 patients in queue at time t
        Z = [Z_0]                 # Number of total patients in ICU at time t
        F = [F_0]                 # Number of infectious patients i ICU at time t
        H = [H_0]                 # Number of active HCP at time t
        J = [J_0]                 # Number of inactive HCP at time t
        q_D = [q_D_0]             # Number of queue departures at time t (untreated)
        q_D_F = [q_D_F_0]         # Number of queue departures who were infectious
        icu_D = [icu_D_0]         # Number of ICU departures at time t (treated)
        N_T = [n_0]               # Current ICU size maximum at time t


        # While the current time is less than the end time
        while t[-1] < T:

          # Check capacity of the unit based on active HCP
          N_T = check_capacity(n_current, H, N_T, r)

          # Find the time that the next event occurs
          # Check all rates
          lmbda_t = check_lambda(t[-1], lmbda)
          mu_t = check_mu(Z[-1],F[-1],N_T[-1], mu_1, mu_2)
          rho_t = check_rho(X_occ, rho_1, rho_2)
          eta_t = check_eta(F[-1],H[-1], eta)
          nu_t = check_nu(J[-1], nu)

          # Total rate
          rate = lmbda_t + mu_t + rho_t + eta_t + nu_t

          # Time update
          delta = np.random.exponential(scale=1/rate,size=1)
          t.append(float(t[-1] + delta))

          # Determine which event occurs at updated time
          events = [1, 2, 3, 4, 5]
          probs = [lmbda_t/rate, rho_t/rate, mu_t/rate, eta_t/rate, nu_t/rate]
          Y = np.random.choice(a=events, p=probs)

          # Event occurence
          if Y == 1:
            X_occ, Z, F, H, J, q_D, q_D_F, icu_D = arrival(t[-1], N_T[-1], X_occ, Z, F, H, J, q_D, q_D_F, icu_D)
          elif Y == 2:
            X_occ, Z, F, H, J, q_D, q_D_F, icu_D = q_depart(X_occ, Z, F, H, J, q_D, q_D_F, icu_D)
          elif Y == 3:
            X_occ, Z, F, H, J, q_D, q_D_F, icu_D = icu_depart(N_T[-1], X_occ, Z, F, H, J, q_D, q_D_F, icu_D)
          elif Y == 4:
            X_occ, Z, F, H, J, q_D, q_D_F, icu_D = HCP_infect(X_occ, Z, F, H, J, q_D, q_D_F, icu_D)
          elif Y == 5:
            X_occ, Z, F, H, J, q_D, q_D_F, icu_D = HCP_return(X_occ, Z, F, H, J, q_D, q_D_F, icu_D)
          else:
            print("Something went wrong. :c")

          # Determine length of queue
          #X.append(len(X_occ))

          # Determine number of infectious patients in queue
          #C.append(sum(X_occ))

        # Record run
        #run_t.append(t)
        #run_X.append(X)
        #run_C.append(C)
        #run_Z.append(Z)
        #run_F.append(F)
        #run_H.append(H)
        run_q_D.append(q_D)
        #run_q_D_F.append(q_D_F)
        #run_N_T.append(N_T)

      # Calculating average accumulation & stdev of untreated deaths for the given pair
      mean, stdev = get_mortality_stats(run_q_D)
      U_mort_avg[n_current,H_current] = mean
      U_mort_std[n_current,H_current] = stdev

  return U_mort_avg, U_mort_std

#------------------------------------------------------------------------------------------------
# Sensitivity analysis...

#-----------------------------------------------------------------------------------------------
# Testing area
single_model(14, 30, lmbda, mu_1, mu_2, 1/10, 1/20, eta, nu, r)
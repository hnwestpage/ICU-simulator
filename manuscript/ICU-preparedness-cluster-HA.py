# ICU PREPAREDNESS SIMULATOR 2.1-------------------------------------

# imports
import numpy as np
import math
import time

#---------------------------------------------------------------------
# PRELIMINARY FUNCTIONS
# Budget operations calculator

# Function to find linear combinations of HCP and beds given a fixed budget
def get_budget_options(U, a, b):
  # Rescale budget based on constant costs
  U = (1-c-d)*U

  # Get max values for budget
  HCP_max = int(U/a)
  n_max = int(U/b)

  budget_array = []
  # Build array of pairs
  for i in range(n_max):
    inner_vec = []
    for j in range(HCP_max):
      pair = (i+1,j+1)
      inner_vec.append(pair)
    budget_array.append(inner_vec)

  return budget_array

# Parameter initializations

# generate random number for saving file
num = math.floor((10**6)*np.random.rand())
num_str = str(num)

# Resources                # Description
U = 2.6 * 10**7            # Hospital budget in USD (total expense = $26 million)
a = 4.0 * 10**5            # Average HCP yearly salary ($400,000)
b = 1.50 * 10**5           # Average cost of one ICU bed for a year ($150,000)
c = 0.32                   # Percent of budget spent on constant costs (Utilities, maintenence)
d = 0.28                   # Percent of budget spent on nursing/tech personnel

# Get combinations based on allotted budget
#budget_array = get_budget_options(U, a, b)

# Rates
lmbda = 9         # Baseline arrival rate of patients to ICU; 9.5 ~ {McManus (2004), Begen (2024)} or 
mu_1 = 1/3.4         # Departure rate of baseline (non infectious) patients from ICU (1/recovery); Moira et al., 2017
mu_2 = 1/8.0         # Departure rate of COVID-19 (infectious) patients from ICU (1/recovery)
rho_1 = 0.12         # Departure rate of baseline patients from the queue (renege)
rho_2 = rho_1        # Departure rate of COVID-19 patients from the queue (renege)
eta = 0.008          # Rate at which HCP becomes infected in the workplace
nu = 1/5           # Rate at which infected HCP recover/return to work

# Proportions
mort_ICU = 0.245           # Proportion of baseline patients who perish in the ICU- Kim et al. (2021), Vincent et al. (2020)
mort_Q = 0.507             # Proportion of baseline patients who perish in the queue- Boumendil et al. (2012)
mort_ICU_19 = 0.245        # Proportion of COVID-19 patients who perish in the ICU- !NEEDS REF
mort_Q_19 = 0.507          # Proportion of COVID-19 patients who perish in the queue- !NEEDS REF
crit = 0.0112              # Proportion of infected indiviuals requiring ICU- Moon et al., (2022)
r =  1/(9.3)               # Threshold ratio of HCP to beds (servers)

# SEIR parameters - pathogen like SARS-CoV-2 B.1.1.529 (Omicron)
beta_0 = 0               # contact rate between susceptibles and infectives (NO EPIDEMIC)
gamma = 1/3              # reciprocal of mean exposed period (3 days)
alpha = 1/5              # reciprocal of mean recovery period (5 days)
tt = 6                   # number of time steps per day (~4 hour time steps)
N = 72046                # population size (High resource setting = 72,046)

# Miscellaneous
M = 200                           # Number of sample paths
T = 365                           # Time (days)
t_ints = list(range(0,T,1))       # Averaging function time intervals
x = np.linspace(0, T, T*tt + 1)   # SEIR timescale variable

# Print out current parameters--------------------------------------------------------------------------------------------------------------------------------
print('Scenario HA')
print('\nJob ID: ',num)
print("\nResources:-----------------\nBudget: ",U,"\nHCP salary: ",a,"\nBed cost: ",b,"\n\% Constant costs: ",c,"\n\% Nursing staff costs: ",d)
print("\nRates:-----------------\nLambda: ",lmbda,"\nmu (baseline): ",mu_1,"\nmu (COVID): ",mu_2,"\nReneging: ",rho_1,"\nHCP infection rate per patient: ",eta,"\nHCP recovery: ",nu)
print("\nProportions:-----------------\n\% Critical condition: ",crit,"\nCoverage ratio (HCP:Patients): ",r)
print("\nSEIR Parameters:-----------------\nbeta: ",beta_0,"\ngamma: ",gamma,"\nalpha: ",alpha,"\nPopulation size: ",N)
print("\nTotal sample paths per scenario: ",M)

# Simulated SEIR-pandemic
# SEIR ODE model normalized for community population
def SEIR_zero_dim(e, i, N):
  # initialize seeds
  
  E_0 = e
  I_0 = i
  R_0 = 0
  S_0 = 1-E_0-I_0

  # Calculating basic reproduction number
  #R0 = beta_0/(alpha*gamma)
  #print(R0)

  # Initialize lists to record number of patients in compartment at each time step
  S = [S_0]
  E = [E_0]
  I = [I_0]
  R = [R_0]

  # Begin for loop to iterate through each time step for duration of model
  for t in range(T*tt):

      # Implementing previous equations
      dS = (-beta_0*S[t]*I[t]) / tt
      dE = (beta_0*S[t]*I[t] - gamma*E[t]) / tt
      dI = (gamma*E[t] - alpha*I[t]) / tt
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
      
      # 1 = S + E + I + R
      S.append(1 - E[t] - I[t] - R[t])

  return I

# Set and Get Functions ###########################################################################################

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


"""# Function to check the current queue departure rate
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

  return rho_1_t + rho_2_t"""

# Function to check the current queue departure rate
def check_rho(X_occ, rho_1):
  i = len(X_occ)

  if i <= 0:
    return 0
  else:
    rho_1_t = i*rho_1

  return rho_1_t


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

# Event Functions ####################################################################################################

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

# Finding times and averages ######################################################################################

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

# function to average the number of queue departures/deaths with patient types
def get_mortality_stats(run_q_D, run_q_D_F):
  entries_1 = []                          # empty list for grabbing last entry of run
  entries_2 = []                          # empty list for grabbing last entry of run

  for i in range(M):                    # for each run
    entries_1.append(run_q_D[i][-1])      # append the last entry of the run; counting ALL queue abandoning events

  for i in range(M):                    # for each run
    entries_2.append(run_q_D_F[i][-1])      # append the last entry of the run; only counting COVID patient queue abandoning events

  # "entries" denotes non-COVID-19 departures, "entries_F" denotes COVID-19 departures
  entries = [a - b for a, b in zip(entries_1, entries_2)]   # counting only NON-COVID ptient queue abandonments
  entries_F = entries_2

  # Get the average of each
  avg_1 = np.mean(entries)    # Non-COVID
  avg_2 = np.mean(entries_F)  # ONLY COVID
  avg_3 = np.mean(entries_1)  # Total

  # Get standard deviation of each
  std_1 = np.std(entries) # non covid
  std_2 = np.std(entries_F) # covid
  std_3 = np.std(entries_1) # total
  
  """ # Commented out, as we are not designating the reason or proportion of queue abandonments due to mortality
  U_mort_mean_1 = mort_Q * avg_1
  U_mort_mean_19 = mort_Q_19 * avg_2

  U_mort_std_1 = mort_Q * std_1
  U_mort_std_19 = mort_Q_19 * std_2

  U_mort_mean = U_mort_mean_1 + U_mort_mean_19
  U_mort_std = U_mort_std_1 + U_mort_std_19
  """
  return avg_1, avg_2, avg_3, std_1, std_2, std_3

# function to average the number of queue departures/deaths SIMPLIFIED
def get_mortality_stats_simp(run_q_D):
  entries_1 = []                          # empty list for grabbing last entry of run

  for i in range(M):                    # for each run
    entries_1.append(run_q_D[i][-1])      # append the last entry of the run; counting ALL queue abandoning events

  # Get the average of each
  avg_3 = np.mean(entries_1)  # Total

  # Get standard deviation of each
  std_3 = np.std(entries_1) # total
  
  """ # Commented out, as we are not designating the reason or proportion of queue abandonments due to mortality
  U_mort_mean_1 = mort_Q * avg_1
  U_mort_mean_19 = mort_Q_19 * avg_2

  U_mort_std_1 = mort_Q * std_1
  U_mort_std_19 = mort_Q_19 * std_2

  U_mort_mean = U_mort_mean_1 + U_mort_mean_19
  U_mort_std = U_mort_std_1 + U_mort_std_19
  """
  return avg_3, std_3

# Function to find the (beds/hcp combo) minimizer of queue deaths
def find_mort_min(U_mort_avg):
  x_star = np.argmin(U_mort_avg)
  return x_star

#---------------------------------------------------------------------
# SIMULATOR LOOP

# SEIR model
# Seeding S, E, I, R
e = 0
i = 0

# Implmenting model: % exposed, % infected, community size
I_normal = SEIR_zero_dim(e, i, N)

# List comprehension method for scaling infection wave by community population size
I = [x * N for x in I_normal]

#**************************************************************************************************************************************
# SINGLE SCENARIO MODEL****************************************************************************************************************
def single_model(H, n, lmbda, mu_1, mu_2, rho_1, eta, nu, ratio):

  # Single scenario cell- fixed HCP and n
  H_current = H
  n_current = n

  # Print baseline traffic density lambda/(n x mu)
  traffic = lmbda/(n_current*mu_1)
  #print('The baseline traffic density is '+str(traffic)+'%')

  # Sensitivity analysis***
  run_q_D = []
  run_q_D_F = []

  #run_rho_t_tracer = []

  for i in range (M):
    run_t = []
    #rho_t_tracer = []

    # print("Run number",str(i+1))

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
      rho_t = check_rho(X_occ, rho_1)
      eta_t = check_eta(F[-1],H[-1], eta)
      nu_t = check_nu(J[-1], nu)

      #rho_t_tracer.append(rho_t)

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

      run_t.append(t)

    # Printing resulting array of runs*** Single Scenario Testing Only
    """print('Timestamp list:',t)
    print('Number in queue',X)
    print('Number of infectious patients in queue',C)
    print('Number in ICU',Z)
    print('Number of infectious patients in ICU',F)
    print('Number of active clinicians',H)
    print('Queue departures',q_D)
    print('Infectious queue departures',q_D_F)
    print('Effective ICU Capacity',N_T)"""

    # Sensitivity analysis***
    run_q_D.append(q_D)
    run_q_D_F.append(q_D_F)

    #run_rho_t_tracer.append(rho_t_tracer)
  
  """#averaging run of rho_t
  avg_rho_t = find_avg_run(run_t,run_rho_t_tracer)
  print(avg_rho_t)"""

  # Return the mortality stats*** Sensitivity Analysis Runs
  mean1, mean2, mean3, stdev1, stdev2, stdev3 = get_mortality_stats(run_q_D, run_q_D_F)

  # non-COVID-19 patients
  U_mort_avg = mean1
  U_mort_std = stdev1

  # COVID-19 patients
  U_mort_avg_19 = mean2
  U_mort_std_19 = stdev2

  # total or ALL patients
  Total_avg = mean3
  Total_stdev = stdev3

  return [float(U_mort_avg), float(U_mort_std), float(U_mort_avg_19), float(U_mort_std_19), float(Total_avg), float(Total_stdev)]

    
#****************************************************************************************************************************************
# MULTI SCENARIO MODEL*******************************************************************************************************************
def multi_model(U, a, b, mu_1, mu_2, eta, nu):
  # track run time
  start_time = time.time()

  # Get combinations based on allotted budget
  budget_array = get_budget_options(U, a, b)

  # Calculate portion of budget for HCP and Beds
  U = (1-c-d)*U

  # Get HCP and bed maxes
  max_HCP = int(U/a)
  max_beds = int(U/b)

  # Printing parameters
  print('Budget = ',U,', HCP salary = ',a,', Bed cost = ',b,', max HCP = ',max_HCP,', max beds = ',max_beds)

  # Initializing logs of untreated deaths per budget option
  U_mort_avg = np.zeros((max_beds,max_HCP))      # dimension (n x H)
  U_mort_std = np.zeros((max_beds,max_HCP))    # dimension (n x H)

  # Per budget scenario
  for row in budget_array:
    for pair in row:
      # Establishing the current scenario
      n_current = pair[0]
      H_current = pair[1]

      # Empty list for each run of q_D
      run_q_D = []

      # Pre-pandemic traffic density
      traffic = lmbda/(n_current*mu_1)

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
          rho_t = check_rho(X_occ, rho_1)
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
        run_q_D.append(q_D)
        #run_q_D_F.append(q_D_F)

        # Printing completed run
        #print('Run ',i,' of ',M,' complete for ',n_current,' beds and ',H_current,' HCP.')

      # Calculating average accumulation & stdev of untreated deaths for the given pair
      mean, stdev = get_mortality_stats_simp(run_q_D)
      U_mort_avg[n_current-1,H_current-1] = mean
      U_mort_std[n_current-1,H_current-1] = stdev

  # Save abandonment counts as an array
  U_mort_avg_array = np.asarray(U_mort_avg)
  U_mort_std_array = np.asarray(U_mort_std)

  # Export arrays as csv files
  np.savetxt("Mean_abandonment_array_HA_"+num_str+".csv",U_mort_avg_array,delimiter=",")
  np.savetxt("Std_dev_abandonment_array_HA_"+num_str+".csv", U_mort_std_array, delimiter=",")

  """print("Mortality Averages")
  for row in U_mort_avg:
    print(row)

  print("Mortality standard deviations")
  for row in U_mort_std:
    print(row)"""

  print("\n--- %s seconds ---" % (time.time() - start_time))


# Running Multi-Model (Total budget, HCP annual salary, annual ICU bed cost, 1/LOS baseline, 1/LOS epidemic, HPC transmisison, HCP recovery)
multi_model(U, a, b, mu_1, mu_2, eta, nu)

"""#-----------------------------------------------------------------------------------------------
# SENSITIVITY ANALYSIS**********************************************************************************************
# Parameter settings and ranges

# Rural/Urban setting [(5,10)/(15,35)]
H = 5
n = 10

# Parameter ranges
lmbda_range = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7]   # Begen et al., (2024) {URBAN} Murthy et al. (2015) {low income/rural}
mu_1_range =  [1/2, 1/3, 1/4, 1/5, 1/6, 1/7, 1/8, 1/9, 1/10, 1/11, 1/12, 1/13, 1/14, 1/15, 1/16, 1/17]  # Moitra et al., (2017)
mu_2_range =  [1/4, 1/5, 1/6, 1/7, 1/8, 1/9, 1/10, 1/11, 1/12, 1/13, 1/14, 1/15] # Essafi et al., 2022
rho_1_range =  [0.1, 0.08, 0.06, 0.04, 0.02, 0.01]            # Lo, (2001)
#rho_2_range = [] # NEEDS REFERENCE (hard to find, absorbing rho2 into rho1)
eta_range = [0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011]  #Lai et al., 2021
nu_range = [1/7, 1/8, 1/9, 1/10]                              #ACEP n.d. (accessed 2025)
ratio_range =  [1/8, 1/9, 1/10, 1/11]                         # Bhatla & Ryskina (2020)

# Resulting abandonment averages (non-covid-19)
lmbda_deaths_avg = []
mu_1_deaths_avg = []
mu_2_deaths_avg = []
rho_1_deaths_avg = []
#rho_2_deaths_avg = []
eta_deaths_avg = []
nu_deaths_avg = []
ratio_deaths_avg = []

# Resulting mortality standard deviations (non-covid-19)
lmbda_deaths_std = []
mu_1_deaths_std = []
mu_2_deaths_std = []
rho_1_deaths_std = []
#rho_2_deaths_std = []
eta_deaths_std = []
nu_deaths_std = []
ratio_deaths_std = []

# Resulting abandonment averages (Covid-19)
lmbda_deaths_avg_F = []
mu_1_deaths_avg_F = []
mu_2_deaths_avg_F = []
rho_1_deaths_avg_F = []
#rho_2_deaths_avg_F = []
eta_deaths_avg_F = []
nu_deaths_avg_F = []
ratio_deaths_avg_F = []

# Resulting mortality standard deviations (Covid-19)
lmbda_deaths_std_F = []
mu_1_deaths_std_F = []
mu_2_deaths_std_F = []
rho_1_deaths_std_F = []
#rho_2_deaths_std = []
eta_deaths_std_F = []
nu_deaths_std_F = []
ratio_deaths_std_F = []

# Resulting abandonment averages TOTALS
lmbda_deaths_avg_T = []
mu_1_deaths_avg_T = []
mu_2_deaths_avg_T = []
rho_1_deaths_avg_T = []
#rho_2_deaths_avg_T = []
eta_deaths_avg_T = []
nu_deaths_avg_T = []
ratio_deaths_avg_T = []

# Resulting mortality standard deviations TOTALS
lmbda_deaths_std_T = []
mu_1_deaths_std_T = []
mu_2_deaths_std_T = []
rho_1_deaths_std_T = []
#rho_2_deaths_std = []
eta_deaths_std_T = []
nu_deaths_std_T = []
ratio_deaths_std_T = []

for each in lmbda_range:
  # Run the model   H, n, lmbda, mu_1, mu_2, rho_1, rho_2, eta, nu, ratio
  death_avg1, death_stdev1, death_avg2, death_stdev2, death_avg_Total, death_stdev_Total = single_model(H, n, each, mu_1, mu_2, rho_1, rho_2, eta, nu, r)

  # Record average deaths in list
  lmbda_deaths_avg.append(death_avg1)
  # Record standard deviation in list
  lmbda_deaths_std.append(death_stdev1)

  # Record average deaths in list
  lmbda_deaths_avg_F.append(death_avg2)
  # Record standard deviation in list
  lmbda_deaths_std_F.append(death_stdev2)

  # Record average deaths in list
  lmbda_deaths_avg_T.append(death_avg_Total)
  # Record standard deviation in list
  lmbda_deaths_std_T.append(death_stdev_Total)


# Sensitivity on mu_1
for each in mu_1_range:
  # Run the model
  death_avg1, death_stdev1, death_avg2, death_stdev2, death_avg_Total, death_stdev_Total = single_model(H, n, lmbda, each,  mu_2, rho_1, rho_2, eta, nu, r)

  # Record average deaths in list
  mu_1_deaths_avg.append(death_avg1)
  # Record standard deviation in list
  mu_1_deaths_std.append(death_stdev1)

  # Record average deaths in list
  mu_1_deaths_avg_F.append(death_avg2)
  # Record standard deviation in list
  mu_1_deaths_std_F.append(death_stdev2)

  # Record average deaths in list
  mu_1_deaths_avg_T.append(death_avg_Total)
  # Record standard deviation in list
  mu_1_deaths_std_T.append(death_stdev_Total)


# Sensitivity on mu_2
for each in mu_2_range:
  # Run the model
  death_avg1, death_stdev1, death_avg2, death_stdev2, death_avg_Total, death_stdev_Total = single_model(H, n, lmbda, mu_1, each, rho_1, rho_2, eta, nu, r)

  # Record average deaths in list
  mu_2_deaths_avg.append(death_avg1)
  # Record standard deviation in list
  mu_2_deaths_std.append(death_stdev1)

  # Record average deaths in list
  mu_2_deaths_avg_F.append(death_avg2)
  # Record standard deviation in list
  mu_2_deaths_std_F.append(death_stdev2)

  # Record average deaths in list
  mu_2_deaths_avg_T.append(death_avg_Total)
  # Record standard deviation in list
  mu_2_deaths_std_T.append(death_stdev_Total)


# Sensitivity on rho_1
for each in rho_1_range:
  # Run the model
  death_avg1, death_stdev1, death_avg2, death_stdev2, death_avg_Total, death_stdev_Total = single_model(H, n, lmbda, mu_1, mu_2, each, rho_2, eta, nu, r)

  # Record average deaths in list
  rho_1_deaths_avg.append(death_avg1)
  # Record standard deviation in list
  rho_1_deaths_std.append(death_stdev1)

  # Record average deaths in list
  rho_1_deaths_avg_F.append(death_avg2)
  # Record standard deviation in list
  rho_1_deaths_std_F.append(death_stdev2)

  # Record average deaths in list
  rho_1_deaths_avg_T.append(death_avg_Total)
  # Record standard deviation in list
  rho_1_deaths_std_T.append(death_stdev_Total)


# Sensitivity on rho_2 COMMENTED OUT
#for each in rho_2_range:
  # Run the model
  #death_avg, death_stdev = single_model(H, n, lmbda, mu_1, mu_2, rho_1, each, eta, nu, r)

  # Record average deaths in list
  #rho_2_deaths_avg.append(death_avg)
  # Record standard deviation in list
  #rho_2_deaths_std.append(death_stdev)


# Sensitivity on eta
for each in eta_range:
  # Run the model
  death_avg1, death_stdev1, death_avg2, death_stdev2, death_avg_Total, death_stdev_Total = single_model(H, n, lmbda, mu_1, mu_2, rho_1, rho_2, each, nu, r)

  # Record average deaths in list
  eta_deaths_avg.append(death_avg1)
  # Record standard deviation in list
  eta_deaths_std.append(death_stdev1)

  # Record average deaths in list
  eta_deaths_avg_F.append(death_avg2)
  # Record standard deviation in list
  eta_deaths_std_F.append(death_stdev2)

  # Record average deaths in list
  eta_deaths_avg_T.append(death_avg_Total)
  # Record standard deviation in list
  eta_deaths_std_T.append(death_stdev_Total)


# Sensitivity on nu
for each in nu_range:
  # Run the model
  death_avg1, death_stdev1, death_avg2, death_stdev2, death_avg_Total, death_stdev_Total = single_model(H, n, lmbda, mu_1, mu_2, rho_1, rho_2, eta, each, r)

  # Record average deaths in list
  nu_deaths_avg.append(death_avg1)
  # Record standard deviation in list
  nu_deaths_std.append(death_stdev1)

  # Record average deaths in list
  nu_deaths_avg_F.append(death_avg2)
  # Record standard deviation in list
  nu_deaths_std_F.append(death_stdev2)

  # Record average deaths in list
  nu_deaths_avg_T.append(death_avg_Total)
  # Record standard deviation in list
  nu_deaths_std_T.append(death_stdev_Total)


# Sensitivity on coverage ratio
for each in ratio_range:
  # Run the model
  death_avg1, death_stdev1, death_avg2, death_stdev2, death_avg_Total, death_stdev_Total = single_model(H, n, lmbda, mu_1, mu_2, rho_1, rho_2, eta, nu, each)

  # Record average deaths in list
  ratio_deaths_avg.append(death_avg1)
  # Record standard deviation in list
  ratio_deaths_std.append(death_stdev1)

  # Record average deaths in list
  ratio_deaths_avg_F.append(death_avg2)
  # Record standard deviation in list
  ratio_deaths_std_F.append(death_stdev2)

  # Record average deaths in list
  ratio_deaths_avg_T.append(death_avg_Total)
  # Record standard deviation in list
  ratio_deaths_std_T.append(death_stdev_Total)

# appending all lists into one
# means
aggregate_list = []
aggregate_list.append(lmbda_deaths_avg)
aggregate_list.append(mu_1_deaths_avg)
aggregate_list.append(mu_2_deaths_avg)
aggregate_list.append(rho_1_deaths_avg)
#aggregate_list.append(rho_2_deaths_avg)
aggregate_list.append(eta_deaths_avg)
aggregate_list.append(nu_deaths_avg)
aggregate_list.append(ratio_deaths_avg)

# standard devs
aggregate_list_std = []
aggregate_list_std.append(lmbda_deaths_std)
aggregate_list_std.append(mu_1_deaths_std)
aggregate_list_std.append(mu_2_deaths_std)
aggregate_list_std.append(rho_1_deaths_std)
#aggregate_list_std.append(rho_2_deaths_std)
aggregate_list_std.append(eta_deaths_std)
aggregate_list_std.append(nu_deaths_std)
aggregate_list_std.append(ratio_deaths_std)

# ----------------------------------------------------------
# means
aggregate_list_F = []
aggregate_list_F.append(lmbda_deaths_avg_F)
aggregate_list_F.append(mu_1_deaths_avg_F)
aggregate_list_F.append(mu_2_deaths_avg_F)
aggregate_list_F.append(rho_1_deaths_avg_F)
#aggregate_list.append(rho_2_deaths_avg)
aggregate_list_F.append(eta_deaths_avg_F)
aggregate_list_F.append(nu_deaths_avg_F)
aggregate_list_F.append(ratio_deaths_avg_F)

# standard devs
aggregate_list_std_F = []
aggregate_list_std_F.append(lmbda_deaths_std_F)
aggregate_list_std_F.append(mu_1_deaths_std_F)
aggregate_list_std_F.append(mu_2_deaths_std_F)
aggregate_list_std_F.append(rho_1_deaths_std_F)
#aggregate_list_std.append(rho_2_deaths_std)
aggregate_list_std_F.append(eta_deaths_std_F)
aggregate_list_std_F.append(nu_deaths_std_F)
aggregate_list_std_F.append(ratio_deaths_std_F)

print("Mortality averages",aggregate_list_F)
print("Mortality standard deviations",aggregate_list_std_F)

# ----------------------------------------------------------
# means
aggregate_list_T = []
aggregate_list_T.append(lmbda_deaths_avg_T)
aggregate_list_T.append(mu_1_deaths_avg_T)
aggregate_list_T.append(mu_2_deaths_avg_T)
aggregate_list_T.append(rho_1_deaths_avg_T)
#aggregate_list.append(rho_2_deaths_avg)
aggregate_list_T.append(eta_deaths_avg_T)
aggregate_list_T.append(nu_deaths_avg_T)
aggregate_list_T.append(ratio_deaths_avg_T)

# standard devs
aggregate_list_std_T = []
aggregate_list_std_T.append(lmbda_deaths_std_T)
aggregate_list_std_T.append(mu_1_deaths_std_T)
aggregate_list_std_T.append(mu_2_deaths_std_T)
aggregate_list_std_T.append(rho_1_deaths_std_T)
#aggregate_list_std.append(rho_2_deaths_std)
aggregate_list_std_T.append(eta_deaths_std_T)
aggregate_list_std_T.append(nu_deaths_std_T)
aggregate_list_std_T.append(ratio_deaths_std_T)

print("Mortality averages",aggregate_list_T)
print("Mortality standard deviations",aggregate_list_std_T)
#-----------------------------------------------------------------------------------------------"""
"""# Running sensitivity analysis

#N_range = [5*(10**4), 7.5*(10**4), 1*(10**5), 1.25*(10**5), 1.5*(10**5)]

ren_list_50k_High = single_model(15, 43, 9, 1/(3.4), 1/10, 1/20, 0.008, 1/8.5, 1/9.5)
ren_list_50k_Low = single_model(5, 17, 3.5, 1/(3.4), 1/10, 1/20, 0.008, 1/8.5, 1/9.5)

# Convert lists into arrays
a1 = np.array(ren_list_50k_High)
a2 = np.array(ren_list_50k_Low)

# Save arrays as CSV files
np.savetxt("Community size 150K High Arrivals.csv",a1,delimiter=",")
np.savetxt("Community size 150K Low Arrivals.csv",a2,delimiter=",")

# Progress updater
print("Completed")
#print("General averages\n",death_avg1,"\n", death_stdev1,"\n", death_avg2,"\n", death_stdev2, "\n", death_avg_Total, "\n", death_stdev_Total)"""
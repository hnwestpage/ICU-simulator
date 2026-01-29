import numpy as np
import math
import matplotlib.pyplot as plt

lmbda = 8
mu = 1/3.3
tau = 1
n = 30
r = lmbda/mu

list_1 = []
for k in range(n):
    num = (np.float_power(r,k))/(math.factorial(k))
    list_1.append(num)

# probability of losing a patient WITHOUT reneging    
p_l = ((np.float_power(r,n))/math.factorial(n))/(sum(list_1))
p_s = 1 - p_l

#print("probability of being treated (no reneging)= ",p_s)

# probability of losing a patient WITH reneging
def find_p_s_tau(lmbda, mu, n):
    r = lmbda/mu

    if r != n:
        p_lr = (((np.float_power(r,n))/math.factorial(n))*np.exp(-mu*tau*(n - r)))/(sum(list_1)+ (((np.float_power(r,n))/math.factorial(n))*(r*np.exp(-mu*tau*(n - r))-n)/(r-n)))
    else:
        p_lr = ((n**n)/(math.factorial(n)))/(sum(list_1)+ (n**n)/(math.factorial(n))*(1+n*mu*tau))

    p_sr = 1 - p_lr
    return p_sr

# Graphing the probability of being served with different arrival rates

lmbda_list = np.linspace(0,15)

p_s_list_30 = []
n = 30
for l in lmbda_list:
    p_s_list_30.append(find_p_s_tau(l, mu, n))

p_s_list_25 = []
n = 25
for l in lmbda_list:
    p_s_list_25.append(find_p_s_tau(l, mu, n))

p_s_list_20 = []
n = 20
for l in lmbda_list:
    p_s_list_20.append(find_p_s_tau(l, mu, n))

p_s_list_15 = []
n = 15
for l in lmbda_list:
    p_s_list_15.append(find_p_s_tau(l, mu, n))

plt.plot(lmbda_list, p_s_list_30, label="30 beds")
plt.plot(lmbda_list, p_s_list_25, label="25 beds")
plt.plot(lmbda_list, p_s_list_20, label="20 beds")
plt.plot(lmbda_list, p_s_list_15, label="15 beds")
plt.xlabel("Number of arrivals per day")
plt.ylabel("Probability of being treated")
plt.legend()
plt.title("Probability of being treated; LOS = 3.3 days")
plt.show()


# Graphing the probability of being served with different disacharge rates

los_list = np.linspace(0.05,21)
mu_list = []
for i in los_list:
    mu_list.append(1/i)

p_s_list_30 = []
n = 30
for m in mu_list:
    p_s_list_30.append(find_p_s_tau(lmbda, m, n))

p_s_list_25 = []
n = 25
for m in mu_list:
    p_s_list_25.append(find_p_s_tau(lmbda, m, n))

p_s_list_20 = []
n = 20
for m in mu_list:
    p_s_list_20.append(find_p_s_tau(lmbda, m, n))

p_s_list_15 = []
n = 15
for m in mu_list:
    p_s_list_15.append(find_p_s_tau(lmbda, m, n))

plt.plot(los_list, p_s_list_30, label="30 beds")
plt.plot(los_list, p_s_list_25, label="25 beds")
plt.plot(los_list, p_s_list_20, label="20 beds")
plt.plot(los_list, p_s_list_15, label="15 beds")
plt.xlabel("Mean length of stay")
plt.ylabel("Probability of being treated")
plt.legend()
plt.title("Probability of being treated; 8 arrivals per day")
plt.show()

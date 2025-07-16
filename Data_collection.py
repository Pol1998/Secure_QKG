##################
#Packages required
##################
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.special import comb
import os

################################
#Mathematical functions required
################################
#Photon loss, probability of "k" losses
def p_loss(k, N, eta):
    """Probability of loss with binomial-like distribution."""
    if k < 0:
        return 0 
    else:
        return comb(N, k) * ((1 - eta) ** k) * (eta ** (N - k))

#Dark counts, probability of "k" dark counts
def p_dark(k,lambd):
    """Dark count probability using Poisson distribution."""
    if k < 0:
        return 0  
    return ((lambd) ** k * np.exp(-lambd)) / math.factorial(k)

def f_coherence(x, gamma):
    return math.exp(-0.5 * gamma * x ** 2)

#Auxilliary function
def var_theta(x, y, z, eta, gamma, lambd):
    if z >= x:
        output = f_coherence(x, gamma) * p_dark(y,lambd) * np.sqrt(p_loss(2*x+y-z, 2*x, eta) * p_loss(2*x+y-z, x, eta))
    else:
        output = f_coherence(z-y, gamma) * p_dark(y,lambd) * np.sqrt(p_loss(2*x+y-z, 2*x, eta) * p_loss(x, x, eta))
    return 1/2 * output

#Objective function
def Phi(N, eta, gamma, lambd, l, s):
    output = 0
    if s == l-N and l >= N:
        suma = 0
        for m in range(s+1):
            suma += var_theta(N, m, l, eta, gamma, lambd)
        output = -suma
    elif l < N:
        output = -var_theta(N, s, l, eta, gamma, lambd)
    return output

###############
#Task functions
###############
def bf_optimizer(N, eta, gamma, lambd):
    u_bound = 1 + int( + 2 * N + N * lambd + 5 * np.sqrt(N * lambd)) #Sets an upper bound within the search space for l and s. (5 standard deviations deviated from the mean dark counted number applied to the maximum initial number of photons, i,e, 2N).
    cost = 10 ** (-5) #Slightly greater than 0 to avoid the unwanted outcome opt_l = opt_s = 0.
    opt_l = 0
    opt_s = 0
    for l in range(u_bound):
        for s in range(u_bound -1):
            if l > s:
                cost_aux = Phi(N, eta, gamma, lambd, l, s)
                if cost_aux < cost:
                    cost = cost_aux
                    opt_l = l
                    opt_s = s
    return [opt_l, opt_s, cost] 

###################################
#Collect data for plot 1 [L,S,W](N)
###################################
#Noiseless implementation
data_noiseless = []
eta = 1 - 10 ** (-5)
gamma = 0
lambd = 0
for N in range(100):
    data_noiseless.append(bf_optimizer(N+1, eta, gamma, lambd))

#Estimated-noise implementation
data_estimated = []
eta = 7*10 ** (-1) #Estimated value for a 50km-long fiber.
gamma = 10 ** (-4) #Arbitrarily set value (because of lack of experimental evidence on what an appropiate value would be)
lambd = 10 ** (-3) #Value acknowledged within our cited source.
for N in range(100):
    data_estimated.append(bf_optimizer(N+1, eta, gamma, lambd))

#Large-noise
data_noisy = []
eta = 10 ** (-5) #Estimated value for a 100km-long fiber.
gamma = 1 #Arbitrarily set value (because of lack of experimental evidence on what an appropiate value would be)
lambd = 10 ** (-1) #Value acknowledged within our cited source.
for N in range(100):
    data_noisy.append(bf_optimizer(N+1, eta, gamma, lambd))

###############################
#Collect data for plot 2 W(eta)
###############################
#Anomalous dependence on the transmitivity by the coherence-witness value W

data_anomalous = [] 
N = 10
gamma = 10 ** (-3) #Arbitrarily set value (because of lack of experimental evidence on what an appropiate value would be)
lambd = 10 ** (-3) #Value acknowledged within our cited source.
for i in range(101):
    data_anomalous.append([i, bf_optimizer(N, i/100 + 10 ** (-5), gamma, lambd)[2]])

#Note: save data in a convenient format
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.special import comb
import os

# === Mathematical functions ===

def p_loss(k, N, eta):
    """Probability of losing k photons out of N with efficiency eta."""
    if k < 0 or k > N:
        return 0
    return comb(N, k) * ((1 - eta) ** k) * (eta ** (N - k))

def p_dark(k, lambd):
    """Poisson probability of k dark counts with mean lambda."""
    if k < 0:
        return 0
    return (lambd ** k * np.exp(-lambd)) / math.factorial(k)

def f_coherence(x, gamma):
    """Coherence decay factor based on parameter gamma."""
    return math.exp(-0.5 * gamma * x ** 2)

def var_theta(x, y, z, eta, gamma, lambd):
    """Helper term used in computing the cost function."""
    if y - z + x < 0:
        return 0
    return p_dark(y, lambd) * np.sqrt(
        p_loss(x + y - z, 2 * x, eta) * p_loss(x + y - z, x, eta)
    )

def Phi(N, eta, gamma, lambd, l):
    """
    Cost function Phi.
    l is a threshold parameter; this function returns how 'bad' it is.
    """
    u_bound = 1 + int(lambd + 5 * np.sqrt(lambd))  # reasonable upper bound
    if l > N + u_bound:
        return 0
    suma = f_coherence(N, gamma) * sum(
        var_theta(N, m, l, eta, gamma, lambd) for m in range(l + 1)
    )
    return -0.5 * suma

def p_counts(N, m, lambd, eta):
    """Total probability of getting m counts (considering both bases)."""
    u_bound = 1 + int(lambd + 5 * np.sqrt(lambd))
    suma_1 = 0
    for i in range(2 * N + 1):
        for j in range(u_bound):
            if 2 * N - i + j == m:
                suma_1 += p_loss(i, 2 * N, eta) * p_dark(j, lambd)
    suma_2 = 0
    for i in range(N + 1):
        for j in range(u_bound):
            if N - i + j == m:
                suma_2 += p_loss(i, N, eta) * p_dark(j, lambd)
    return 0.5 * (suma_1 + suma_2)

def p_conditional(m, N, f, lambd, eta):
    """Conditional probability of count m given bit f (0 or 1)."""
    u_bound = 1 + int(lambd + 5 * np.sqrt(lambd))
    suma_1 = sum(
        p_loss(i, 2 * N, eta) * p_dark(j, lambd)
        for i in range(2 * N + 1)
        for j in range(u_bound)
        if 2 * N - i + j == m
    )
    suma_2 = sum(
        p_loss(i, N, eta) * p_dark(j, lambd)
        for i in range(N + 1)
        for j in range(u_bound)
        if N - i + j == m
    )
    if suma_1 == 0 and suma_2 == 0:
        return 0.5  # avoid division by zero
    return suma_1 / (suma_1 + suma_2) if f == 0 else suma_2 / (suma_1 + suma_2)

def QBER(N, lambd, eta):
    """Quantum Bit Error Rate (QBER) estimation."""
    u_bound = 1 + int(lambd + 5 * np.sqrt(lambd))
    error_0, error_1 = 0, 0
    for m in range(2 * N + u_bound):
        pc0 = p_conditional(m, N, 0, lambd, eta)
        pc1 = p_conditional(m, N, 1, lambd, eta)
        if pc1 > pc0:
            error_0 += p_counts(N, m, lambd, eta) * pc0
        else:
            error_1 += p_counts(N, m, lambd, eta) * pc1
    return error_0 + error_1

def bf_optimizer(N, eta, gamma, lambd):
    """Brute-force optimization: finds best threshold l minimizing Phi."""
    best_cost = 1e-5
    opt_l = 0
    for l in range(N + 1):
        cost_aux = Phi(N, eta, gamma, lambd, l)
        if cost_aux < best_cost:
            best_cost = cost_aux
            opt_l = l
    return [opt_l, best_cost]

# === Output Directory ===
# IMPORTANT: Change this path if you want to save output files elsewhere
# For example: output_dir = "/home/yourname/results"
output_dir = "./output"
os.makedirs(output_dir, exist_ok=True)

# === Generate and Save Plot 1 Data (L vs N) ===

def generate_plot1_data():
    """
    Generates optimal threshold L vs photon number N for different noise settings.
    Saves text files with two columns: N and L.
    """
    settings = {
        "noiseless": (1 - 1e-5, 0, 0),
        "low_noise": (0.9, 1e-4, 1e-3),
        "estimated": (0.7, 1e-4, 1e-3),
        "high_noise": (0.3, 1e-4, 1e-3),
        "noisy": (0.01, 1, 1e-3)
    }

    N = np.arange(1, 101)

    for label, (eta, gamma, lambd) in settings.items():
        data_L = [bf_optimizer(n, eta, gamma, lambd)[0] for n in N]
        np.savetxt(
            os.path.join(output_dir, f"data_{label}.txt"),
            np.column_stack((N, data_L)),
            header="N L", comments=''
        )

# === Generate and Save Plot 2 Data (Cost & QBER vs N) ===

def generate_plot2_data():
    """
    Generates cost function values and QBER vs N for selected noise settings.
    Saves two text files per setting: one for cost, one for QBER.
    """
    settings = {
        "low_noise": (0.9, 1e-4, 1e-3),
        "estimated": (0.7, 1e-4, 1e-3),
        "high_noise": (0.3, 1e-4, 1e-3),
    }

    N = np.arange(1, 101)

    for label, (eta, gamma, lambd) in settings.items():
        data_cost = []
        data_qber = []
        for n in N:
            l_opt, cost = bf_optimizer(n, eta, gamma, lambd)
            data_cost.append(cost)
            data_qber.append(QBER(n, lambd, eta))
        np.savetxt(
            os.path.join(output_dir, f"data_cost_{label}.txt"),
            np.column_stack((N, data_cost)),
            header="N Xi", comments=''
        )
        np.savetxt(
            os.path.join(output_dir, f"data_qber_{label}.txt"),
            np.column_stack((N, data_qber)),
            header="N QBER", comments=''
        )

# === Run Simulations ===

generate_plot1_data()
generate_plot2_data()

print("All data saved successfully to:", output_dir)


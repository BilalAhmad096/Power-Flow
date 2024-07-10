import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import inv

# Define system parameters
n_buses = 33
line_data = [
    (1, 2, 0.0922 + 0.0470j), (2, 3, 0.4930 + 0.2511j), (3, 4, 0.3660 + 0.1864j), (4, 5, 0.3811 + 0.1941j),
    (5, 6, 0.8190 + 0.7070j), (6, 7, 0.1872 + 0.6188j), (7, 8, 0.7114 + 0.2351j), (8, 9, 1.0300 + 0.7400j),
    (9, 10, 1.0440 + 0.7400j), (10, 11, 0.1966 + 0.0650j), (11, 12, 0.3744 + 0.1238j), (12, 13, 1.4680 + 1.1550j),
    (13, 14, 0.5416 + 0.7129j), (14, 15, 0.5910 + 0.5260j), (15, 16, 0.7463 + 0.5450j), (16, 17, 1.2890 + 1.7210j),
    (17, 18, 0.7320 + 0.5740j), (2, 19, 0.1640 + 0.1565j), (19, 20, 1.5042 + 1.3554j), (20, 21, 0.4095 + 0.4784j),
    (21, 22, 0.7089 + 0.9373j), (3, 23, 0.4512 + 0.3083j), (23, 24, 0.8980 + 0.7091j), (24, 25, 0.8960 + 0.7011j),
    (6, 26, 0.2030 + 0.1034j), (26, 27, 0.2842 + 0.1447j), (27, 28, 1.0590 + 0.9337j), (28, 29, 0.8042 + 0.7006j),
    (29, 30, 0.5075 + 0.2585j), (30, 31, 0.9744 + 0.9630j), (31, 32, 0.3105 + 0.3619j), (32, 33, 0.3410 + 0.5302j)
]  # (Bus from, Bus to, Line Impedance (ohm))

# Load data
loads = [
    (1, 0.000, 0.000), (2, 0.100, 0.060), (3, 0.090, 0.040), (4, 0.120, 0.080), (5, 0.060, 0.030), 
    (6, 0.060, 0.020), (7, 0.200, 0.100), (8, 0.200, 0.100), (9, 0.060, 0.020), (10, 0.060, 0.020), 
    (11, 0.045, 0.030), (12, 0.060, 0.035), (13, 0.060, 0.035), (14, 0.120, 0.080), (15, 0.060, 0.010), 
    (16, 0.060, 0.020), (17, 0.060, 0.020), (18, 0.090, 0.040), (19, 0.090, 0.040), (20, 0.090, 0.040), 
    (21, 0.090, 0.040), (22, 0.090, 0.040), (23, 0.090, 0.050), (24, 0.420, 0.200), (25, 0.420, 0.200), 
    (26, 0.060, 0.025), (27, 0.060, 0.025), (28, 0.060, 0.020), (29, 0.120, 0.070), (30, 0.200, 0.600), 
    (31, 0.150, 0.070), (32, 0.210, 0.100), (33, 0.060, 0.040)
] # (Bus No, Active Power Demand(MW), Reactive Power Demand(MVAR))

# Distributed generation data
dist_gen = [(1, 4.0, 2.5), (18, 0.2, 0), (22, 0.2, 0), (25, 0.2, 0), (33, 0.2, 0)]  # (Bus No, Active Power Capacity (MW), Reactive Power Capacity(MVAR))

# Power demands
P_loads = np.array([load[1] * 1000 for load in loads])  # Convert to kW
Q_loads = np.array([load[2] * 1000 for load in loads])  # Convert to kVAR

# Slack bus voltage
Vslack = 253  # Slack voltage in volts
Deltaslack = 0 # Slack bus angle in degrees

# Convert line data to Y-bus matrix
Y_bus = np.zeros((n_buses, n_buses), dtype=complex)
for (i, j, z) in line_data:
    Y_bus[i-1, i-1] += 1/z
    Y_bus[j-1, j-1] += 1/z
    Y_bus[i-1, j-1] -= 1/z
    Y_bus[j-1, i-1] -= 1/z

# Initialize bus voltages and angles
V = np.ones(n_buses) * Vslack
Delta = np.zeros(n_buses)

# Initialize power mismatch
P_spec = P_loads / 1000  # Convert kW to MW
Q_spec = Q_loads / 1000  # Convert kVAR to MVAR

# Add distributed generation
for bus, P_gen, Q_gen in dist_gen:
    P_spec[bus-1] -= P_gen
    Q_spec[bus-1] -= Q_gen

# Tolerance and max iterations
tolerance = 1e-6
max_iter = 10

for iteration in range(max_iter):
    P_calc = np.zeros(n_buses)
    Q_calc = np.zeros(n_buses)
    
    for i in range(n_buses):
        for j in range(n_buses):
            P_calc[i] += V[i] * V[j] * (Y_bus[i, j].real * np.cos(Delta[i] - Delta[j]) + Y_bus[i, j].imag * np.sin(Delta[i] - Delta[j]))
            Q_calc[i] += V[i] * V[j] * (Y_bus[i, j].real * np.sin(Delta[i] - Delta[j]) - Y_bus[i, j].imag * np.cos(Delta[i] - Delta[j]))

    # Power mismatch
    dP = P_spec[1:] - P_calc[1:]
    dQ = Q_spec[1:] - Q_calc[1:]
    mismatch = np.r_[dP, dQ]
    
    # Check convergence
    if np.max(np.abs(mismatch)) < tolerance:
        break
    
    # Form the Jacobian matrix
    J = np.zeros((2*(n_buses-1), 2*(n_buses-1)))
    for i in range(1, n_buses):
        for j in range(1, n_buses):
            if i == j:
                for k in range(n_buses):
                    if k != i:
                        J[i-1, j-1] += V[i] * V[k] * (-Y_bus[i, k].real * np.sin(Delta[i] - Delta[k]) + Y_bus[i, k].imag * np.cos(Delta[i] - Delta[k]))
                J[i-1, j-1] -= V[i]**2 * Y_bus[i, i].imag
                J[i-1, n_buses-1+j-1] = V[i] * Y_bus[i, i].real + V[i] * np.sum([V[k] * (Y_bus[i, k].real * np.cos(Delta[i] - Delta[k]) + Y_bus[i, k].imag * np.sin(Delta[i] - Delta[k])) for k in range(n_buses) if k != i])
            else:
                J[i-1, j-1] = V[i] * V[j] * (Y_bus[i, j].real * np.sin(Delta[i] - Delta[j]) - Y_bus[i, j].imag * np.cos(Delta[i] - Delta[j]))
                J[i-1, n_buses-1+j-1] = V[i] * (Y_bus[i, j].real * np.cos(Delta[i] - Delta[j]) + Y_bus[i, j].imag * np.sin(Delta[i] - Delta[j]))
    
    for i in range(1, n_buses):
        for j in range(1, n_buses):
            if i == j:
                for k in range(n_buses):
                    if k != i:
                        J[n_buses-1+i-1, j-1] += V[i] * V[k] * (Y_bus[i, k].real * np.cos(Delta[i] - Delta[k]) + Y_bus[i, k].imag * np.sin(Delta[i] - Delta[k]))
                J[n_buses-1+i-1, j-1] -= V[i]**2 * Y_bus[i, i].real
                J[n_buses-1+i-1, n_buses-1+j-1] = V[i] * Y_bus[i, i].imag + V[i] * np.sum([V[k] * (-Y_bus[i, k].real * np.sin(Delta[i] - Delta[k]) + Y_bus[i, k].imag * np.cos(Delta[i] - Delta[k])) for k in range(n_buses) if k != i])
            else:
                J[n_buses-1+i-1, j-1] = V[i] * V[j] * (-Y_bus[i, j].real * np.cos(Delta[i] - Delta[j]) - Y_bus[i, j].imag * np.sin(Delta[i] - Delta[j]))
                J[n_buses-1+i-1, n_buses-1+j-1] = V[i] * (-Y_bus[i, j].real * np.sin(Delta[i] - Delta[j]) + Y_bus[i, j].imag * np.cos(Delta[i] - Delta[j]))
    
    # Solve for voltage and angle corrections
    dx = np.linalg.solve(J, mismatch)
    
    # Update voltage angles and magnitudes
    Delta[1:] += dx[:n_buses-1]
    V[1:] += dx[n_buses-1:]
    
# print("Voltages (V):", V)
# print("Angles (degrees):", np.degrees(Delta))

# Print voltages and angles in the specified format
print("Bus_No    Voltage (V)    Angle (degrees)")
for i in range(len(V)):
    print(f"{i+1:<10} {V[i]:<15.8f} {np.degrees(Delta[i]):<18.8f}")

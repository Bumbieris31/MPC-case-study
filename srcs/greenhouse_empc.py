import numpy as np
from casadi import *
import scipy.io
import pandas as pd
from scipy.interpolate import interp1d


def vanHentenODEs(x, u, d, p):
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28 = p

    w1 = 0; w2 = 0; w3 = 0; w4 = 0  # No disturbances

    x1 = x[0]; x2 = x[1]; x3 = x[2]; x4 = x[3]

    u1 = u[0]  # Supply rate of CO2 in mg*m^{-2}*s^{-1}
    u2 = u[1]  # Ventilation rate through the vents in mm s^{-1}
    u3 = u[2]  # Energy supply by the heating system in W*m^{-2}

    d1 = d[0]  # Incoming radiation in W m^{-2}
    d2 = d[1]  # Outside CO2 in kg m^{-3}
    d3 = d[2]  # Outdoor temperature in ^oC
    d4 = d[3]  # Outdoor humidity content C_{H2O} in kg m^{-3}

    phi = p4*d1 + (-p5*x3**2 + p6*x3 - p7)*(x2 - p8)
    PhiPhot_c = (1-np.exp(-p3*x1))*(p4*d1*(-p5*x3**2 + p6*x3 - p7)*(x2-p8))/phi  # gross canopy photosynthesis rate
    PhiVent_c = (u2*10**(-3) + p11)*(x2-d2)  # mass exchange of CO2 through the vents
    PhiVent_h = (u2*10**(-3) + p11)*(x4 - d4)  # canopy transpiration
    PhiTransp_h = p21*(1 - np.exp(-p3*x1))*(p22/(p23*(x3+p24))*np.exp(p25*x3/(x3+p26))-x4)  # mass exchange of H2O through the vents

    dx1dt = (p1*PhiPhot_c - p2*x1*2**(x3/10 - 5/2))*(1+w1)
    dx2dt = 1/p9*(-PhiPhot_c + p10*x1*2**(x3/10 - 5/2) + u1*10**(-6) - PhiVent_c)*(1+w2)
    dx3dt = 1/p16*(u3 - (p17*u2*10**(-3) + p18)*(x3 - d3) + p19*d1)*(1+w3)
    dx4dt = 1/p20*(PhiTransp_h - PhiVent_h)*(1+w4)

    dxdt = np.array([dx1dt, dx2dt, dx3dt, dx4dt])
    return dxdt


# Function to perform explicit rk4
def f_rk4(x, u, d, Delta, p):
    n_rk4 = 1
    delta_rk4 = Delta / n_rk4
    for i in range(n_rk4):
        k_1 = vanHentenODEs(x, u, d, p)
        k_2 = vanHentenODEs(x + 0.5 * delta_rk4 * k_1, u, d, p)
        k_3 = vanHentenODEs(x + 0.5 * delta_rk4 * k_2, u, d, p)
        k_4 = vanHentenODEs(x + delta_rk4 * k_3, u, d, p)
        x = x + (1 / 6) * (k_1 + 2 * k_2 + 2 * k_3 + k_4) * delta_rk4
    return x


# Function to convert states to outputs
def h_meas(x, p):
    y_1 = 10 ** (3) * x[0]  # Weight in g m^{-2}
    y_2 = p[11] * (x[2] + p[12]) / (p[13] * p[14]) * x[1] * 10 ** (3)  # Indoor CO2 in ppm 10^{3}
    y_3 = x[2]  # Air Temp in  ^oC
    y_4 = p[11] * (x[2] + p[12]) / (11 * np.exp(p[26] * x[2] / (x[2] + p[27]))) * x[3] * 10 ** (2)  # RH C_{H2O} in %
    return np.array([y_1, y_2, y_3, y_4])


# load parameters (p)
p = np.ravel(scipy.io.loadmat('vanHentenParams.mat')["p"])
# load weather (weather)
weather = scipy.io.loadmat('weatherConditions.mat')["weather"][0]
# data = weather["weather"]


# time step
Delta = 3 * 300  # 15 minutes
N = 12

# initial conditions
t0 = 0
x0 = np.array([0.0035, 0.001, 15, 0.008])

# Downsample weather from 5 mins to 15 mins
raw_data = np.array([weather['iGlob'], weather['co2'], weather['tOut'], weather['hum']])
d_cl = np.array([np.ravel(raw_data[0][0]), np.ravel(raw_data[1][0]), np.ravel(raw_data[2][0]), np.ravel(raw_data[3][0])])

# d_cl = np.interp(d_cl, (1, (Delta / 300), len(d_cl)))
x_weather = np.arange(1, d_cl.shape[1] + 1, Delta/300)
d_cl_interp = np.zeros((d_cl.shape[0], len(x_weather)))

for i in range(d_cl.shape[0]):
    d_cl_interp[i, :] = np.interp(x_weather, np.arange(d_cl.shape[1]), d_cl[i, :])



# Define cost parameters
c_co2 = 0.42 * Delta  # per kg/s of CO2
c_q = 6.35E-9 * Delta  # per W of heat
c_dw = -16  # price per kg of dry weight

# Define optimization problem
opti = casadi.Opti()

# State and input trajectory
xhat = opti.variable(4, N + 1)
uhat = opti.variable(3, N)
# Initial state
x0_p = opti.parameter(4)
# Weather in trajectory
d_p = opti.parameter(4, N)

cost = 0
opti.subject_to(xhat[:, 0] == x0_p)
for k in range(N):
    # Enforce model
    opti.subject_to(xhat[:, k + 1] == f_rk4(xhat[:, k], uhat[:, k], d_p[:, k], Delta, p))

    # Input constraints
    opti.subject_to(0 <= uhat[0, k] <= 1.2)
    opti.subject_to(0 <= uhat[1, k] <= 7.5)
    opti.subject_to(0 <= uhat[2, k] <= 150)

    # State/output constraints
    yhat = h_meas(xhat[:, k], p)
    opti.subject_to(0 <= yhat[0])
    opti.subject_to(0 <= yhat[1] <= 1.6)
    opti.subject_to(10 <= yhat[2] <= 25)
    opti.subject_to(0 <= yhat[3] <= 80)

    # Add stage cost
    cost = cost + c_dw * (xhat[0, k + 1] - xhat[0, k]) + c_co2 * uhat[0, k] * 10 ** (-6) + c_q * uhat[2, k]

# Define objective function
opti.minimize(cost)

# Solver options
solver_options = {}
# Suppress output
solver_options['ipopt.print_level'] = 0
solver_options['print_time'] = 0
solver_options['ipopt.sb'] = 'yes'
# Select solver
opti.solver('ipopt', solver_options)

# Closed-loop simulation
Nsim = 4 * 24  # Number of simulation time steps

# Initialize closed-loop state and input variables
x_cl = np.zeros((4, Nsim + 1))
y_cl = np.zeros((4, Nsim + 1))
u_cl = np.zeros((3, Nsim))

# Set initial condition
x_cl[:, 0] = x0
y_cl[:, 0] = h_meas(x0, p)
print("Simulation step:")
total_cost = 0
for k in range(Nsim):
    print(k)
    # Set current state at time k as initial state in MPC optimization
    opti.set_value(x0_p, x_cl[:, k])
    opti.set_value(d_p, d_cl[:, t0 + k:t0 + k + N])
    # Solve optimization problem
    sol = opti.solve()
    # Save warm-start
    opti.set_initial(sol.value_variables())
    # Get optimal input trajectory
    u_opt = sol.value(uhat)
    # Take only the *first* input
    u_cl[:, k] = u_opt[:, 0]
    # Simulate system
    x_cl[:, k + 1] = f_rk4(x_cl[:, k], u_cl[:, k], d_cl[:, t0 + k], Delta, p)
    y_cl[:, k + 1] = h_meas(x_cl[:, k + 1], p)

    total_cost = total_cost + c_dw * (x_cl[0, k + 1] - x_cl[0, k]) + c_co2 * u_cl[0, k] * 10 ** (-6) + c_q * u_cl[2, k]

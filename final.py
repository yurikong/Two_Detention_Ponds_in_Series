import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import math

plt.rc('text', usetex=True)         # use latex
plt.rc('font', family='serif')      # set fonts


# inflow function q_i,1(t): water inflow discharge rate
# convert unit of time from hours to seconds to match the unit of q
def q_i_1(t):
    if 3600 <= t < 12600:               # from t0 to t1
        return 1 / 4500 * (t - 3600)    # left line with positive slope
    elif 12600 <= t < 21600:            # from t1 to t2
        return -1 / 4500 * (t - 21600)  # right line with negative slope
    else:                               # 0 otherwise
        return 0


# system of ODE
# params (unit)
# D1 = D2 = outlet pipe diameters (m)
# A1 = bottom area of pond 1 (m^2)
# A2 = bottom area of pond 2 (m^2)
# H1 = max depth of pond 1 (m)
# H2 = max depth of pond 2 (m)
# h1_0 = h2_0 = initial water depth in both ponds (0)
def f(y, t, D1, D2, A1, A2, H1, H2):
    h1 = y[0]
    h2 = y[1]
    if h1 > H1:     # when water in pond 1 exceeds its max height
        h1 = H1     # excessive water will be lost (cap at max height of pond 1)
    if h2 > H2:     # when water in pond 2 exceeds its max height
        h2 = H2     # excessive water will be lost (cap at max height of pond 2)
    if h1 < 0:      # physically, height of water in pond 1 cannot be less than 0
        h1 = 0      # less than empty = empty
    if h2 < 0:      # physically, height of water in pond 2 cannot be less than 0
        h2 = 0      # less than empty = empty
    dh1_dt = (q_i_1(t) - np.pi * D1 ** 2 / 4 * math.sqrt(2 * 9.8 * h1)) / A1
    dh2_dt = (np.pi * D1 ** 2 / 4 * math.sqrt(2 * 9.8 * h1) - np.pi * D2 ** 2 / 4 * math.sqrt(2 * 9.8 * h2)) / A2
    return [dh1_dt, dh2_dt]


# outflow function q_o,2(t): water outflow discharge rate
def q_o_2(h2, D2):
    return np.pi * D2 ** 2 / 4 * math.sqrt(2 * 9.8 * h2)


# values os params
D1 = D2 = 0.2
A1 = 2000
A2 = 1000
H1 = 5
H2 = 4
h1_0 = h2_0 = 0
# time spans from 0 to 80 hours
# use seconds instead of hours to calculate
t = np.linspace(0, 80 * 3600 + 1, 1000)
IC = [h1_0, h2_0]       # initial ponds depths
H = integrate.odeint(f, IC, t, args=(D1, D2, A1, A2, H1, H2))   # solve the system of ODE
h1, h2 = H.T
# excess water in each pond will be capped at the max height of the pond
# physically, water level cannot be lower than 0
for i in range(len(t)):
    if h1[i] < 0:
        h1[i] = 0
    if h2[i] < 0:
        h2[i] = 0
    if h1[i] > H1:
        h1[i] = H1
    if h2[i] > H2:
        h2[i] = H2
# plotting first graph: depth vs time
plt.plot(t / 3600, h1, 'r-', label='$h_1(t)$')      # use hours instead of seconds in x-axis
plt.plot(t / 3600, h2, 'b--', label='$h_2(t)$')     # use hours instead of seconds in x-axis
plt.title('Pond Water Depth vs. Time')
plt.xlabel('Time (hr)')
plt.ylabel('Depth (m)')
plt.legend()        # show the legend
plt.show()
plt.close()
# plotting second graph: water discharge rates
plt.plot(t / 3600, [q_i_1(i) for i in t], 'r-', label='$q_{i,1}(t)$')           # use hours instead of seconds in x-axis
plt.plot(t / 3600, [q_o_2(i, D2) for i in h2], 'b--', label='$q_{o,2}(t)$')     # use hours instead of seconds in x-axis
plt.title('Pond Water Discharge Rates')
plt.xlabel('Time (hr)')
plt.ylabel('Discharge Rate ($m^3/s$)')
plt.legend()
plt.show()

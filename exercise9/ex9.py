# ======================================================================
# Theoretical and Computational Star Formation
# Assignment sheet 9
# Clara Kretzschmar
# Timon Danowski
# Tom Carron
# ======================================================================
import numpy as np

# definfing computational domain and shock initial conditions (given q1, q2, shock location and number of particles N)
def initial(a, b, N, q_1, q_2, shock_loc):
    q_array = np.zeros((N, 5))
    x = np.zeros(N)
    d_x = dx(a, b, N)
    shock_loc = N / b * shock_loc
    # initializing the grid
    for i in range(N):
        x[i] = d_x * i
        if i < shock_loc:
            # initial conditions
            for j in range(5):
                q_array[i][j] = q_1[j]
        else:
            for j in range(5):
                q_array[i][j] = q_2[j]
    return q_array


# Delta x, Domain a,b. N particles.
def dx(a, b, N):
    dx = (b - a) / N
    return dx


# function to determine the energy from omega and gamma, from the ideal gas equation
def E(w, gamma):
    v = 0
    for i in range(1, 3):
        v += w[i] ** 2
    return w[4] / (gamma - 1) + 0.5 * w[0] * v


def max_eigenvalue(w, gamma=1.4):
    """
    compute the maximum absolute eigenvalue of the system
    input:  primitive variables (density, velocity_x, velocity_y, velocity_z, pressure) and gamma
    see exercise sheet 7
    """
    #If statements to avoid possible division by zero
    if w[0] != 0:
        lambda_1 = w[1]  # v_x
        lambda_2 = w[1] - np.sqrt(
            gamma * w[4] / w[0]
        )  # v_x - sqrt(gamma * pressure / density)
        lambda_3 = w[1] + np.sqrt(
            gamma * w[4] / w[0]
        )  # v_x + sqrt(gamma * pressure / density)
    else:
        lambda_1 = w[1]  # v_x
        lambda_2 = w[1]  # v_x
        lambda_3 = w[1]  # v_x


    lambda_max = np.max([np.abs(lambda_1), np.abs(lambda_2), np.abs(lambda_3)])

    return lambda_max


# converting omega in q and f
# w = [rho, vx, vy, vz, p]
def q_function(w, gamma):
    rho = w[0]
    q = np.zeros(5)
    q[0] = rho
    for i in range(1, 3):
        q[i] = rho * w[i]
    q[4] = E(w, gamma)
    return q


# to get w from q
def invert_q(q, gamma):
    f = np.zeros(5)
    if q[0] != 0:
        f[0] = q[0]
        f[1] = q[1] / q[0]
        f[2] = q[2] / q[0]
        f[3] = q[3] / q[0]
        f[4] = (gamma - 1) * (q[4] - 0.5 * q[0] * (f[1] ** 2 + f[2] ** 2 + f[3] ** 2))
    else:
        for i in range(5):
            f[i]=0
    
    return f

#To calculate f_x array from primitive variables
def fx(w, gamma):
    rho = w[0]
    v = 0
    for i in range(1, 3):
        v += w[i] ** 2
    fx = np.zeros(5)
    fx[0] = rho * w[1]
    fx[1] = fx[0] * w[1] + w[4]
    fx[2] = fx[1] * w[2]
    fx[3] = fx[1] * w[3]
    fx[4] = w[1] * (0.5 * rho * np.abs(v) + (gamma * w[4]) / (gamma - 1))
    return fx


# finding the max (absolute value) of two eigenvalues
def lambda_max(lambda_i, lambda_mi):
    return np.maximum(np.amax(lambda_i), np.amax(lambda_mi))


# function to determine the next timestep input lambda l, with the CFL condition: C = 0.001; l = lambda * dt/dx
def getTimestep(l, dx):
    max_step=0.001
    C = 0.5
    cfl=(dx * C / l)
    if cfl < max_step:
        dt=cfl
    else:
        dt=max_step
    return dt


# function to calculate the f_array
def getRiemannFlux(q_array, gamma):
    f = np.zeros((len(q_array), 5))
    for i in range(len(q_array)):
        f[i] = fx(invert_q(q_array[i], gamma), gamma)
    return f


# function to update the Solution
def update(q_array, finterface, dt, dx, l):
    # boundary condition
    new_q_array = np.zeros_like(q_array)
    for j in range(len(q_array)-1):
        if j == 0:
            i = j + 1
        elif j == len(q_array)-1:
            i = j - 1
        else:
            i = j
        temp1 = approx_Riemann_solver(
            finterface[i + 1], finterface[i], q_array[i + 1], q_array[i], l[i + 1], l[i]
        )
        temp2 = approx_Riemann_solver(
            finterface[i], finterface[i - 1], q_array[i], q_array[i - 1], l[i], l[i - 1]
        )
        new_q_array[j] = q_array[j] - (dt / dx * (temp1 - temp2))
    return new_q_array


# calculating fx_i-1/2
# postion i: i, position i-1: im, lambda l
def approx_Riemann_solver(fx_i, fx_im, q_i, q_im, l_i, l_im):
    temp = np.zeros(5)
    l_max = lambda_max(l_i, l_im)
    for i in range(5):
        temp[i] = 0.5 * (fx_im[i] + fx_i[i]) - 0.5 * l_max * (q_i[i] - q_im[i])
    return temp

def convert_var(array_prim_var, gamma=1.4):
    ''' 
    converts primitive variables (density, velocity_x, velocity_y, velocity_z, pressure)
    to conservative variables (density, density * v_x, density * v_y, density * v_x, energy)
    see exercise sheet 8, eq. (9)
    '''
    array_conserv_var = np.zeros(5)
    array_conserv_var[0] = array_prim_var[0]
    array_conserv_var[1] = array_prim_var[0] * array_prim_var[1]
    array_conserv_var[2] = array_prim_var[0] * array_prim_var[2]
    array_conserv_var[3] = array_prim_var[0] * array_prim_var[3]
    v_sq = array_prim_var[1]**2 + array_prim_var[2]**2 + array_prim_var[3]**2 # velocity squared
    array_conserv_var[4] = 0.5 * array_prim_var[0] * v_sq + array_prim_var[4] / (gamma-1.0)

    return array_conserv_var

"""
Function to run the simultaion given boundaries a,b,shock location, number of particles N, adiabatic index gamma,
primitive variables on either side of the shock w1 and w2, and simulation end time tmax. This is the only function that needs to be called to run the sim.
(except for calling initial() to plot the initial shock conditions).
"""
def run(N, a, b, shock_loc, gamma, w_1, w_2, tmax):
    t = 0
    q = initial(a, b, N, convert_var(w_1, gamma), convert_var(w_2, gamma), shock_loc)
    # calculate dx
    d_x = dx(a, b, N)
    counter = 0
    while t < tmax:
        counter += 1
        # calculate the maximum time step
        l = []
        for i in range(len(q)):
            l = np.append(l, max_eigenvalue(invert_q(q[i],gamma), gamma))
            dt = getTimestep(np.max(l), d_x)
        # calculate the RiemannFlux
        finterface = getRiemannFlux(q, gamma)
        # update the q'
        q_new = update(q, finterface, dt, d_x, l)
        q = q_new
        t += dt
        # print(dt)
    return q

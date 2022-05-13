# ======================================================================
# Theoretical and Computational Star Formation
# Assignment sheet 4
# Tom Carron, Clara Kretzschmar, Timon Danowski
# ======================================================================
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

parameters = {"font.family": "STIXGeneral", "mathtext.fontset": "stix"}
plt.rcParams.update(parameters)

# ======================================================================

# Read file
def getElemfromFile(filename, N_t_01, N_t_001, N_particles):
    position = np.zeros((N_t_01, N_particles))
    density = np.zeros((N_t_01, N_particles))
    velocity = np.zeros((N_t_01, N_particles))
    totalEnergy = np.zeros((N_t_01, N_particles))
    time = np.zeros(N_t_001)
    pressure = np.zeros((N_t_01, N_particles))
    # position001 = np.zeros((N_t_001, N_particles))
    # density001 = np.zeros((N_t_001, N_particles))
    # velocity001 = np.zeros((N_t_001, N_particles))
    # pressure001 = np.zeros((N_t_001, N_particles))

    with open(filename, "r") as file:
        lines = file.readlines()
        for t, line in enumerate(lines[1:]):
            for p in range(N_particles):
                position[t, p] = line.split()[p + 1]
                density[t, p] = line.split()[p + N_particles + 1]
                velocity[t, p] = line.split()[p + 2 * N_particles + 1]
                # for p in range(N_particles):
                totalEnergy[t, p] = line.split()[p + 3 * N_particles + 1]
                pressure[t, p] = line.split()[p + 4 * N_particles + 1]
                time[t] = line.split()[0]
            # position001[t, p] = line.split()[p]
            # density001[t, p] = line.split()[p+N_particles+1]
            # velocity001[t, p] = line.split()[p+2*N_particles+1]
            # pressure001[t,p] = line.split()[p+4*N_particles+1]

    return (
        position,
        density,
        velocity,
        totalEnergy,
        time,
        pressure,
    )  # , position001, density001, velocity001, pressure001


# ======================================================================

files = [
    "file_0.010000200.txt",
    "file_0.100000200.txt",
    "file_1.000000200.txt",
    "file_0.0100001000.txt",
    "file_0.1000001000.txt",
    "file_1.0000001000.txt",
    "file_0.0100002000.txt",
    "file_0.1000002000.txt",
    "file_1.0000002000.txt",
]
Nums = [200, 200, 200, 1000, 1000, 1000, 2000, 2000, 2000]
alph = [0.01, 0.1, 1.0, 0.01, 0.1, 1.0, 0.01, 0.1, 1.0]
# path_to_files="results/"

# Time of interest
t = 0.2
positions = []
densities = []
velocities = []
energies = []
times = []
pressures = []
for i in range(len(files)):
    current_file = files[i]
    print(current_file)
    N_particles = Nums[i]
    (
        position_1,
        density_1,
        velocity_1,
        totalEnergy_1,
        time,
        pressure_1,
    ) = getElemfromFile(current_file, 1, 1, N_particles)
    positions.append(position_1)
    densities.append(density_1)
    velocities.append(velocity_1)
    energies.append(totalEnergy_1)
    times.append(time)
    pressures.append(pressure_1)
    # df = pd.read_table(current_file)#, names=["time","position","velocity","internal energy"])#,columns=["time","position","velocity","internal energy"])
    # interest=df.loc[df["time"] == str(t)]
    # print(df)

fig, ax = plt.subplots(3, len(files), sharex=False)
fig.set_size_inches(27, 9)


for k in range(len(files)):
    fig, ax = plt.subplots(3, 1, sharex=True)
    ax[0].set_title(r"$\alpha=$" + str(alph[k]) + r", $N_{SPH}=$" + str(Nums[k]))
    ax[0].plot(
        positions[k], densities[k], "o", marker="+", color="blue", label=r"density"
    )
    # ax[0].set_xlabel("position")
    ax[0].set_ylabel("density")
    ax[1].plot(
        positions[k], velocities[k], "o", marker="+", color="blue", label=r"velocity"
    )
    ax[1].set_ylabel("velocity")
    ax[2].plot(
        positions[k], pressures[k], "o", marker="+", color="blue", label=r"pressure"
    )
    ax[2].set_ylabel("pressure")
    ax[2].set_xlabel("position")
    ax[1].set_ylim(-10, 10)
    plt.xlim(0, 1)
    fig.tight_layout()
    plt.savefig("plots/file" + str(k) + ".png", dpi=400)

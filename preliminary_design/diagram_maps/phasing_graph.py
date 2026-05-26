import matplotlib.pyplot as plt

orbit_rev = range(1, 47)
time = [282, 90, 56, 42, 35, 30, 28, 26, 25, 24, 24, 24, 24, 24, 24, 25, 25, 26, 27, 27, 28, 29, 29, 30, 31, 32, 32, 33, 34, 35, 36, 37, 38, 38, 39, 40, 41, 42, 
        43, 44, 45, 46, 47, 47, 48, 49]

plt.plot(orbit_rev, time)
plt.xlabel("Revolutions in phasing orbit")
plt.ylabel("Time [days]")
plt.title("Phasing optimisation")
plt.show()
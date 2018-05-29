import h5py
import numpy as np
import matplotlib.pyplot as plt

fluxes = dict()
average_fluxes = dict()
identifiers = dict()
n_time_steps_average = 100000

with h5py.File('flux_data_new.hdf5', 'r') as f:
    for phantom in f:
        print(phantom)
        fluxes[phantom] = f['/' + phantom + '/n_accum'][:]

        if not average_fluxes:
            time = f['/' + phantom + '/time'][:]
            time_intervals_inv = np.zeros(len(time))
            for i, t in enumerate(time):
                if i < n_time_steps_average:
                    time_intervals_inv[i] = t
                else:
                    time_intervals_inv[i] = t - time[i - n_time_steps_average]

            time_intervals_inv = 1 / time_intervals_inv

        average_fluxes[phantom] = np.array([
            (f - fluxes[phantom][max(0, i - n_time_steps_average)]) * t
            for i, (t, f) in enumerate(zip(time_intervals_inv, fluxes[phantom]))])

    total_output = np.zeros(len(time))

    for phantom in f:
        if not phantom == '2':
            total_output += average_fluxes[phantom]

#TODO: provisional
identifiers['2'] = 'Input'
identifiers['3'] = '1-Divot'
identifiers['4'] = '3-Intermediate'
identifiers['5'] = '2-Near divot'
identifiers['6'] = '4-Gauge'
#

print('Averaging over', 1 / time_intervals_inv[-1], 'seconds.')
f = plt.figure(1)
plt.xlabel('time (s)')
plt.ylabel('accumulated number of particles')
plt.ticklabel_format(style='sci', axis='y', scilimits=(2,0))

for phantom in sorted(average_fluxes):
    plt.plot(time, average_fluxes[phantom], label=identifiers[phantom])

plt.plot(time, total_output, label='output', linestyle='dashed')
ax = plt.gca()
legend = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
          ncol=3, frameon=True)
legend.get_frame().set_edgecolor('black')
f.savefig("accumulated.pdf", bbox_inches='tight')


f = plt.figure(2)
plt.xlabel('time (s)')
plt.ylabel('number of particles per second')
plt.ticklabel_format(style='sci', axis='y', scilimits=(2,0))

for phantom in sorted(average_fluxes):
    plt.plot(time, average_fluxes[phantom], label=identifiers[phantom])

plt.plot(time, total_output, label='output', linestyle='dashed')
expected_flux = np.array([-150000 for x in time])
plt.plot(time,expected_flux,'r', label='expected (steady state)', linestyle='dashed')

ax = plt.gca()
legend = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
          ncol=3, frameon=True)
legend.get_frame().set_edgecolor('black')

f.savefig("averaged.pdf", bbox_inches='tight')

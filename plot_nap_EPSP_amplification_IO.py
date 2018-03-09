from plot_results import *


mpl.rcParams['font.size'] = 14.

output_file_path = 'data/08112017_nap_amplification_IO_DC_soma_stim_trunk.hdf5'

rec_dict = {}
baseline = 10.  # ms
v_th = -52.

with h5py.File(output_file_path, 'r') as f:
    trial = f.itervalues().next()
    if 'dt' in trial.attrs:
        dt = trial.attrs['dt']
    else:
        dt = 0.02
    if 'duration' in trial.attrs:
        duration = trial.attrs['duration']
    else:
        duration = 450.
    if 'equilibrate' in trial.attrs:
        equilibrate = trial.attrs['equilibrate']
    else:
        equilibrate = 250.
    t = np.arange(0., duration, dt)
    start = int((equilibrate - baseline) / dt)
    # end = start + int(100. / dt)
    end = int((duration - 100.)/ dt)
    offset_t = t[start:end] - equilibrate
    for trial in f.itervalues():
        vrest = round(trial.attrs['vrest'])
        if vrest not in rec_dict:
            rec_dict[vrest] = {}
        ttx = trial.attrs['ttx']
        group = 'ttx%i' % int(ttx)
        if group not in rec_dict[vrest]:
            rec_dict[vrest][group] = {}
            rec_dict[vrest][group]['num_syns'] = [0]
            rec_dict[vrest][group]['peak_vm'] = [trial.attrs['vrest']]
            rec_dict[vrest][group]['peak_g_AMPA'] = [0.]
        rec_dict[vrest][group]['num_syns'].append(trial.attrs['num_syns'])
        rec_dict[vrest][group]['peak_g_AMPA'].append(None)
        for rec in trial['rec'].itervalues():
            description = rec.attrs['description']
            this_rec = np.interp(t, trial['time'], rec)[start:end]
            if description == 'soma':
                rec_dict[vrest][group]['peak_vm'].append(min(v_th, np.max(this_rec)))
            elif 'g_AMPA' in description:
                if rec_dict[vrest][group]['peak_g_AMPA'][-1] is None:
                    rec_dict[vrest][group]['peak_g_AMPA'][-1] = this_rec
                else:
                    rec_dict[vrest][group]['peak_g_AMPA'][-1] = np.add(rec_dict[vrest][group]['peak_g_AMPA'][-1],
                                                                       this_rec)
        for i, g_AMPA in enumerate(rec_dict[vrest][group]['peak_g_AMPA']):
            rec_dict[vrest][group]['peak_g_AMPA'][i] = np.max(g_AMPA)

fig, axes = plt.subplots(2, 2, sharex=True, sharey=True)
axes[1][0].set_xlabel('Number of activated synapses')
axes[1][0].set_xlim(0, 40)
for i, (ttx, title) in enumerate(zip(['ttx0', 'ttx1'], ['Control', 'TTX'])):
    axes[i][0].axhline(y=v_th, color='grey', linestyle='--', linewidth=1)
    axes[i][0].set_title(title, fontsize=mpl.rcParams['font.size'])
    axes[i][0].set_ylabel('Peak Vm (mV)')
    axes[i][0].set_ylim(-70., -40.)
    for c, vrest, label in zip(['k', 'r'], sorted(rec_dict.keys()), ['Resting Vm', '5 mV depolarized']):
        indexes = range(len(rec_dict[vrest][ttx]['num_syns']))
        indexes.sort(key=rec_dict[vrest][ttx]['num_syns'].__getitem__)
        num_syns = np.array(rec_dict[vrest][ttx]['num_syns'])[indexes]
        peak_vm = np.array(rec_dict[vrest][ttx]['peak_vm'])[indexes]
        axes[i][0].plot(num_syns, peak_vm, c=c, label=label)
axes[1][0].set_aspect('auto')
axes[0][0].legend(loc='best', frameon=False, framealpha=0.5)
clean_axes(axes)
fig.tight_layout()

fig, axes = plt.subplots(2, 2, sharex=True, sharey=True)
axes[1][0].set_xlabel('Synaptic AMPA-R conductance (nS)')
axes[1][0].set_xlim(0, 25)
for i, (ttx, title) in enumerate(zip(['ttx0', 'ttx1'], ['Control', 'TTX'])):
    axes[i][0].axhline(y=v_th, color='grey', linestyle='--', linewidth=1)
    axes[i][0].set_title(title, fontsize=mpl.rcParams['font.size'])
    axes[i][0].set_ylabel('Peak Vm (mV)')
    axes[i][0].set_ylim(-70., -40.)
    for c, vrest, label in zip(['k', 'r'], sorted(rec_dict.keys()), ['Resting Vm', '5 mV depolarized']):
        indexes = range(len(rec_dict[vrest][ttx]['peak_g_AMPA']))
        indexes.sort(key=rec_dict[vrest][ttx]['peak_g_AMPA'].__getitem__)
        peak_g_AMPA = np.array(rec_dict[vrest][ttx]['peak_g_AMPA'])[indexes]
        peak_vm = np.array(rec_dict[vrest][ttx]['peak_vm'])[indexes]
        axes[i][0].plot(peak_g_AMPA * 1000., peak_vm, c=c, label=label)
        # print title, label
        # print peak_g_AMPA * 1000
        # print peak_vm
axes[1][0].set_aspect('auto')
axes[0][0].legend(loc='best', frameon=False, framealpha=0.5)
clean_axes(axes)
fig.tight_layout()
plt.show()
plt.close()

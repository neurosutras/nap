from plot_results import *


output_filepath = 'data/20180220_nap_EPSC_amplification_soma.hdf5'
# output_filepath = 'data/20180222_nap_EPSC_amplification_soma.hdf5'

rec_dict = {}

baseline = 10.  # ms

with h5py.File(output_filepath, 'r') as f:
    trial = f.itervalues().next()
    if 'dt' in trial.attrs:
        dt = trial.attrs['dt']
    else:
        dt = 0.01
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
        this_vrest = trial.attrs['vrest']
        this_TTX = trial.attrs['ttx']
        this_bar = trial.attrs['bar']
        this_zd = trial.attrs['zd']
        if this_bar:
            if this_TTX:
                group = 'Bar + TTX'
            else:
                group = 'Bar'
        elif this_zd:
            if this_TTX:
                group = 'ZD + TTX'
            else:
                group = 'ZD'
        elif this_TTX:
            group = 'TTX'
        else:
            group = 'Control'
        #if group == 'Bar' and np.isclose(this_vrest, -56., atol=0.2):
        #    continue
        if group not in rec_dict:
            rec_dict[group] = {}
        if 'vrest' not in rec_dict[group]:
            rec_dict[group]['vrest'] = []
        rec_dict[group]['vrest'].append(this_vrest)
        for rec in trial['rec'].itervalues():
            description = rec.attrs['description']
            if description not in rec_dict[group]:
                rec_dict[group][description] = []
            this_rec = np.interp(t, trial['time'], rec)[start:end]
            if description in ['soma', 'dend', 'distal_spine']:
                offset = np.mean(this_rec[int((baseline - 3.) / dt):int((baseline - 1.) / dt)])
                this_rec -= offset
            rec_dict[group][description].append(this_rec)

rec_types = set()
for group in rec_dict.itervalues():
    rec_types.update(group.keys())
rec_types.remove('vrest')

ymin = {}
ymax = {}
for rec_type in rec_types:
    if rec_type in ['soma', 'dend', 'distal_spine']:
        category = 'vm'
    elif rec_type in ['soma_ina', 'dend_ina', 'ais_ina']:
        category = 'ina'
    elif rec_type in ['soma_ika']:
        category = 'ika'
    elif rec_type in ['soma_ih']:
        category = 'ih'
    else:
        category = None
    if category is not None:
        if category not in ymin:
            ymin[category] = 0.
            ymax[category] = 0.
        for group in rec_dict:
            if rec_type in rec_dict[group]:
                for rec in rec_dict[group][rec_type]:
                    ymin[category] = min(ymin[category], np.min(rec))
                    ymax[category] = max(ymax[category], np.max(rec))
for category in ymax:
    if category == 'vm':
        ymax[category] = min(5., ymax[category])
        ymin[category] = -0.1 * ymax[category]
    elif category == 'ika':
        ymax[category] = 0.001
        ymin[category] = -0.1 * ymax[category]
    elif category == 'ina':
        ymin[category] = -0.18
        ymax[category] = 0.1 * abs(ymin[category])
    else:
        ymax[category] = 0.1 * abs(ymin[category])
        ymin[category] *= 1.1


from matplotlib import cm
this_cm = cm.get_cmap()
num_items = len(rec_dict.itervalues().next()['vrest'])
colors = [this_cm(1.*i/float(num_items-1)) for i in range(num_items)]
groups = ['Control', 'Bar', 'ZD', 'TTX', 'Bar + TTX', 'ZD + TTX']

sort_indexes = {}
for rec_type in rec_types:
    fig, axes = plt.subplots(2, 3, sharey=True, sharex=True)
    for i, group in enumerate([group for group in groups if group in rec_dict]):
        if group not in sort_indexes:
            sort_indexes[group] = range(len(rec_dict[group]['vrest']))
            sort_indexes[group].sort(key=rec_dict[group]['vrest'].__getitem__)
            sort_indexes[group].reverse()
            rec_dict[group]['vrest'] = map(rec_dict[group]['vrest'].__getitem__, sort_indexes[group])
        if rec_type in rec_dict[group]:
            row = i / 3
            col = i % 3
            axes[row][col].set_title(group)
            rec_dict[group][rec_type] = map(rec_dict[group][rec_type].__getitem__, sort_indexes[group])
            for j, vrest in enumerate(rec_dict[group]['vrest']):
                if group == 'Control':
                    axes[row][col].plot(offset_t, rec_dict[group][rec_type][j], color=colors[j],
                                        label='%.1f mV' % vrest)
                else:
                    axes[row][col].plot(offset_t, rec_dict[group][rec_type][j], color=colors[j])
    if rec_type in ['soma', 'dend', 'distal_spine']:
        category = 'vm'
        axes[0][0].set_ylabel('EPSP amp (mV)')
        plt.ylim(ymin[category], ymax[category])
        plt.suptitle(rec_type+'_Vm')
    elif rec_type in ['soma_ina', 'dend_ina', 'ais_ina']:
        category = 'ina'
        axes[0][0].set_ylabel('Current (mA/cm2)')
        plt.ylim(ymin[category], ymax[category])
        plt.suptitle(rec_type)
    elif rec_type in ['soma_ika']:
        category = 'ika'
        axes[0][0].set_ylabel('Current (mA/cm2)')
        plt.ylim(ymin[category], ymax[category])
        plt.suptitle(rec_type)
    elif rec_type in ['soma_ih']:
        category = 'ih'
        axes[0][0].set_ylabel('Current (mA/cm2)')
        plt.ylim(ymin[category], ymax[category])
        plt.suptitle(rec_type)
    else:
        axes[0][0].set_ylabel('Current (nA)')
        plt.suptitle(rec_type)
    axes[0][0].legend(loc='best', frameon=False, framealpha=0.5)
    axes[1][0].set_xlabel('Time (ms)')
    clean_axes(axes)
    fig.tight_layout()
    fig.subplots_adjust(top=0.85, hspace=0.5, wspace=0.2)
# plt.show()
# plt.close()


mpl.rcParams['font.size'] = 14.

rec_type = 'soma'

EPSP_amplitude = {}
for group in groups:
    EPSP_amplitude[group] = []
    for trace in rec_dict[group][rec_type]:
        if group == 'Bar' and trace is rec_dict[group][rec_type][0]:
            EPSP_amplitude[group].append(np.max(trace[:int(50./dt)]))
        else:
            EPSP_amplitude[group].append(np.max(trace))

amplification_ratio = {}
hyper_vm = -62.
depol_vm = -56.
for group in EPSP_amplitude:
    values = rec_dict[group]['vrest']
    hyper = min(range(len(values)), key=lambda i: abs(values[i]-hyper_vm))
    depol = min(range(len(values)), key=lambda i: abs(values[i] - depol_vm))
    amplification_ratio[group] = EPSP_amplitude[group][depol] / EPSP_amplitude[group][hyper]

bar_width = 0.8
ordered_groups = groups
ordered_xlabels = ['Control', 'Bar', 'ZD', 'Control', 'Bar', 'ZD']
bar_indexes = np.array([0., 1., 2., 3.15, 4.15, 5.15])
fig, axes = plt.subplots(2, 2, figsize=(8, 7))
group = 'Control'
vals1 = [amplification_ratio[group] for group in ordered_groups]
rects1 = axes[0][0].bar(bar_indexes, vals1, bar_width, color='k')
axes[0][0].set_xticks(bar_indexes+0.2)
axes[0][0].set_xticklabels(ordered_xlabels, rotation=45, ha='right', fontsize=12)
axes[0][0].tick_params(axis='x', width=0)
axes[0][0].set_ylim(0., 3.)
axes[0][0].set_ylabel('Amplification ratio\n(%i mV/%i mV)' % (depol_vm, hyper_vm))
axes[0][0].set_title('Full model', fontsize=mpl.rcParams['font.size'])
axes[0][0].axhline(y=1., color='r', linestyle='--', linewidth=1)
axes[0][0].set_aspect('auto')
clean_axes(axes)
fig.tight_layout()
plt.show()
plt.close()
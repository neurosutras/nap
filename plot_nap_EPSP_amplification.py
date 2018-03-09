from plot_results import *

mpl.rcParams['font.size'] = 14.

output_file_paths = {'soma': 'data/07162017_nap_amplification_DC_soma_stim_trunk.hdf5',
                     'dend': 'data/07162017_nap_amplification_DC_dend_stim_trunk.hdf5'}

rec_dict = {}
syn_counter = {}
raw_amplitude = {}
amplification_ratio = {}
drug_control_ratio = {}

baseline = 10.  # ms

for section, file_path in output_file_paths.iteritems():
    rec_dict[section] = {}
    syn_counter[section] = {}
    with h5py.File(file_path, 'r') as f:    
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
        end = int((duration - 100.)/ dt)
        offset_t = t[start:end] - equilibrate
        for trial in f.itervalues():
            this_vrest = trial.attrs['vrest']
            this_ap5 = trial.attrs['ap5']
            this_ttx = trial.attrs['ttx']
            if this_ap5:
                if this_ttx:
                    group = 'AP5 + TTX'
                else:
                    group = 'AP5'
            else:
                if this_ttx:
                    group = 'TTX'
                else:
                    group = 'Control'
            if group not in rec_dict[section]:
                rec_dict[section][group] = {}
                rec_dict[section][group]['vrest'] = []
                rec_dict[section][group]['recs'] = {}
                syn_counter[section][group] = []
            rec_dict[section][group]['vrest'].append(this_vrest)
            syn_counter[section][group].append(0)
            for rec in trial['rec'].itervalues():
                description = rec.attrs['description']
                if 'spine' in description:
                    description = 'spine'
                    syn_counter[section][group][-1] += 1
                elif 'i_NMDA' in description:
                    description = 'i_NMDA'
                if description not in rec_dict[section][group]['recs']:
                    rec_dict[section][group]['recs'][description] = []
                if description in ['spine', 'i_NMDA'] and len(rec_dict[section][group]['recs'][description]) < \
                        len(rec_dict[section][group]['vrest']):
                    rec_dict[section][group]['recs'][description].append(None)
                this_rec = np.interp(t, trial['time'], rec)[start:end]
                if description in ['soma', 'dend', 'spine']:
                    offset = np.mean(this_rec[int((baseline - 3.) / dt):int((baseline - 1.) / dt)])
                    this_rec -= offset
                if description in ['spine', 'i_NMDA']:
                    if rec_dict[section][group]['recs'][description][-1] is None:
                        rec_dict[section][group]['recs'][description][-1] = this_rec
                    else:
                        rec_dict[section][group]['recs'][description][-1] = \
                            np.add(rec_dict[section][group]['recs'][description][-1], this_rec)
                else:
                    rec_dict[section][group]['recs'][description].append(this_rec)
        for group in rec_dict[section]:
            for description in [description for description in ['spine', 'i_NMDA'] if
                                description in rec_dict[section][group]['recs']]:
                for i in xrange(len(rec_dict[section][group]['recs'][description])):
                    rec_dict[section][group]['recs'][description][i] /= syn_counter[section][group][i]

for section in rec_dict:
    raw_amplitude[section] = {}
    for group in ['Control', 'TTX', 'AP5']:
        raw_amplitude[section][group] = {}
        for description in rec_dict[section][group]['recs']:
            raw_amplitude[section][group][description] = []
            for trace in rec_dict[section][group]['recs'][description]:
                raw_amplitude[section][group][description].append(np.max(np.abs(trace)))

hyper_vm = -62.
depol_vm = -56.
for section in raw_amplitude:
    amplification_ratio[section] = {}
    for group in raw_amplitude[section]:
        amplification_ratio[section][group] = {}
        values = rec_dict[section][group]['vrest']
        #sorted_indexes.sort(key=rec_dict[section][group]['vrest'].__getitem__)
        hyper = min(range(len(values)), key=lambda i: abs(values[i]-hyper_vm))
        depol = min(range(len(values)), key=lambda i: abs(values[i] - depol_vm))
        for description in raw_amplitude[section][group]:
            amplification_ratio[section][group][description] = raw_amplitude[section][group][description][depol] / \
                                                               raw_amplitude[section][group][description][hyper]

for section in amplification_ratio:
    drug_control_ratio[section] = {}
    for group in ['TTX', 'AP5']:
        drug_control_ratio[section][group] = {}
        for description in amplification_ratio[section][group]:
            drug_control_ratio[section][group][description] = amplification_ratio[section][group][description] / \
                                                              amplification_ratio[section]['Control'][description]


"""
bar_width = 0.8
ordered_descriptions = ['soma', 'dend', 'spine', 'ais_ina', 'dend_ina', 'i_NMDA']
ordered_xlabels = ['Soma EPSP', 'Dend EPSP', 'Spine EPSP', 'AIS INa', 'Dend INa', 'INMDA']
bar_indexes = np.arange(len(ordered_descriptions))
fig, axes = plt.subplots(2, 2, figsize=(8, 7))
section = 'soma'
group = 'Control'
vals1 = [amplification_ratio[section][group][description] for description in ordered_descriptions]
rects1 = axes[0][0].bar(bar_indexes, vals1, bar_width, color='k')  # , bottom=0.1)
axes[0][0].set_xticks(bar_indexes+0.4)
axes[0][0].set_xticklabels(ordered_xlabels, rotation=45, ha='right', fontsize=12)
axes[0][0].tick_params(axis='x', width=0)
axes[0][0].set_yscale('log')
#axes[0][0].set_ylim(0.5, 100.)
axes[0][0].set_ylim(0.1, 10.)
axes[0][0].set_ylabel('Amplification ratio\n(%i mV/%i mV)' % (depol_vm, hyper_vm))
# axes[0][0].set_yticks([1, 10, 100])
axes[0][0].set_yticks([1, 10])
#axes[0][0].set_yticklabels([1, 10, 100])
axes[0][0].set_yticklabels([1, 10])
axes[0][0].set_title('Adjusted Vm: Soma', fontsize=mpl.rcParams['font.size'])
axes[0][0].axhline(y=1., color='grey', linestyle='--', linewidth=1)
axes[0][0].set_aspect('auto')
section = 'dend'
vals2 = [amplification_ratio[section][group][description] for description in ordered_descriptions]
rects2 = axes[0][1].bar(bar_indexes, vals2, bar_width, color='k')  # , bottom=0.1)
axes[0][1].set_xticks(bar_indexes+0.4)
axes[0][1].set_xticklabels(ordered_xlabels, rotation=45, ha='right', fontsize=12)
axes[0][1].tick_params(axis='x', width=0)
axes[0][1].set_ylabel('Amplification ratio\n(%i mV/%i mV)' % (depol_vm, hyper_vm))
axes[0][1].set_yscale('log')
#axes[0][1].set_ylim(0.5, 100.)
axes[0][1].set_ylim(0.1, 10.)
#axes[0][1].set_yticks([1, 10, 100])
#axes[0][1].set_yticklabels([1, 10, 100])
axes[0][1].set_yticks([1, 10])
axes[0][1].set_yticklabels([1, 10])
axes[0][1].set_title('Adjusted Vm: Dend', fontsize=mpl.rcParams['font.size'])
axes[0][1].axhline(y=1., color='grey', linestyle='--', linewidth=1)
axes[0][1].set_aspect('auto')

ordered_descriptions = ['soma', 'dend', 'spine', 'soma', 'dend', 'spine', 'i_NMDA']
ordered_groups = ['AP5'] * 3 + ['TTX'] * 4
ordered_xlabels = ['Soma EPSP', 'Dend EPSP', 'Spine EPSP', 'Soma EPSP', 'Dend EPSP', 'Spine EPSP', 'INMDA']
bar_indexes2 = np.arange(len(ordered_descriptions))
section = 'soma'
vals3 = [drug_control_ratio[section][group][description] for group, description in
        zip(ordered_groups, ordered_descriptions)]
rects3 = axes[1][0].bar(bar_indexes2, vals3, bar_width, color='g', edgecolor='k')
for i in xrange(3):
    rects3[i].set_color('r')
    rects3[i].set_edgecolor('k')
axes[1][0].set_xticks(bar_indexes2+0.4)
axes[1][0].set_xticklabels(ordered_xlabels, rotation=45, ha='right', fontsize=12)
axes[1][0].tick_params(axis='x', width=0)
axes[1][0].set_ylim(0., 1.5)
axes[1][0].set_ylabel('Normalized\namplification ratio\n(Drug/Control)')
axes[1][0].set_title('Adjusted Vm: Soma', fontsize=mpl.rcParams['font.size'])
axes[1][0].axhline(y=1., color='grey', linestyle='--', linewidth=1)
axes[1][0].set_aspect('auto')
section = 'dend'
vals4 = [drug_control_ratio[section][group][description] for group, description in
        zip(ordered_groups, ordered_descriptions)]
rects4 = axes[1][1].bar(bar_indexes2, vals4, bar_width, color='g', edgecolor='k')
for i in xrange(3):
    rects4[i].set_color('r')
    rects4[i].set_edgecolor('k')
axes[1][1].set_xticks(bar_indexes2+0.4)
axes[1][1].set_xticklabels(ordered_xlabels, rotation=45, ha='right', fontsize=12)
axes[1][1].tick_params(axis='x', width=0)
axes[1][1].set_ylabel('Normalized\namplification ratio\n(Drug/Control)')
axes[1][1].set_ylim(0., 1.5)
axes[1][1].set_title('Adjusted Vm: Dend', fontsize=mpl.rcParams['font.size'])
axes[1][1].axhline(y=1., color='grey', linestyle='--', linewidth=1)
axes[1][1].set_aspect('auto')
clean_axes(axes)
fig.tight_layout()
plt.show()
plt.close()
"""

bar_width = 0.8
ordered_descriptions = ['soma', 'dend', 'spine', 'soma', 'dend', 'spine']
ordered_sections = ['soma', 'soma', 'soma', 'dend', 'dend', 'dend']
ordered_xlabels = ['Soma', 'Apical dendrite', 'Dendritic spine', 'Soma', 'Apical dendrite', 'Dendritic spine']
bar_indexes = np.array([0., 1., 2., 3.15, 4.15, 5.15])
fig, axes = plt.subplots(2, 2, figsize=(8, 7))
group = 'Control'
vals1 = [amplification_ratio[section][group][description] for description, section in zip(ordered_descriptions,
                                                                                          ordered_sections)]
rects1 = axes[0][0].bar(bar_indexes, vals1, bar_width, color='k')
axes[0][0].set_xticks(bar_indexes+0.2)
axes[0][0].set_xticklabels(ordered_xlabels, rotation=45, ha='right', fontsize=12)
axes[0][0].tick_params(axis='x', width=0)
axes[0][0].set_ylim(0., 3.)
axes[0][0].set_ylabel('Amplification ratio\n(%i mV/%i mV)' % (depol_vm, hyper_vm))
axes[0][0].set_title('Full model', fontsize=mpl.rcParams['font.size'])
axes[0][0].axhline(y=1., color='r', linestyle='--', linewidth=1)
axes[0][0].set_aspect('auto')

ordered_descriptions = ['soma', 'dend', 'spine', 'soma', 'dend', 'spine']
ordered_groups = ['AP5'] * 3 + ['TTX'] * 3
section = 'soma'
vals3 = [drug_control_ratio[section][group][description] for group, description in
        zip(ordered_groups, ordered_descriptions)]
rects3 = axes[1][0].bar(bar_indexes, vals3, bar_width, color='g', edgecolor='k')
for i in xrange(3):
    rects3[i].set_color('r')
    rects3[i].set_edgecolor('k')
axes[1][0].set_xticks(bar_indexes+0.2)
axes[1][0].set_xticklabels(ordered_xlabels, rotation=45, ha='right', fontsize=12)
axes[1][0].tick_params(axis='x', width=0)
axes[1][0].set_ylim(0., 1.5)
axes[1][0].set_ylabel('Normalized amplification ratio\n(Drug/Control)')
axes[1][0].set_title('Full model', fontsize=mpl.rcParams['font.size'])
axes[1][0].axhline(y=1., color='grey', linestyle='--', linewidth=1)
axes[1][0].set_aspect('auto')
clean_axes(axes)
fig.tight_layout()
plt.show()
plt.close()

"""

descriptions = set()
for section in rec_dict:
    for group in rec_dict[section]:
        descriptions.update(rec_dict[section][group]['recs'].keys())

for section in rec_dict:
    for group in ['Control', 'TTX', 'AP5']:
        for description in descriptions:
            if not description in raw_amplitude[section][group]:
                raw_amplitude[section][group][description] = np.zeros_like(rec_dict[section][group]['vrest'])

ordered_descriptions = ['soma', 'dend', 'spine', 'ais_ina', 'dend_ina', 'i_NMDA']
ordered_subtitles = ['Soma EPSP', 'Dend EPSP', 'Spine EPSP', 'AIS INa', 'Dend INa', 'INMDA']
ordered_ylabels = ['EPSP amplitude (mV)', 'EPSP amplitude (mV)', 'EPSP amplitude (mV)', 'Current density (uA/cm2)',
                   'Current density (uA/cm2)', 'Current amplitude (pA)']
for section, suptitle in zip(['soma', 'dend'], ['Soma', 'Dend']):
    fig, axes = plt.subplots(2, 3, figsize=(9, 6), sharex=True)
    fig.suptitle('Adjusted Vm: %s' % suptitle, fontsize=mpl.rcParams['font.size'])
    for i, description in enumerate(ordered_descriptions):
        row = i / 3
        col = i % 3
        ymax = 0.
        for c, group, label in zip(['k', 'r', 'g'], ['Control', 'AP5', 'TTX'], ['Control', 'AP5', 'Low TTX']):
            indexes = range(len(rec_dict[section][group]['vrest']))
            indexes.sort(key=rec_dict[section][group]['vrest'].__getitem__)
            this_vrest = np.array(rec_dict[section][group]['vrest'])[indexes]
            this_vals = np.array(raw_amplitude[section][group][description])[indexes]
            if description in ['ais_ina', 'dend_ina', 'i_NMDA']:
                this_vals *= 1000.
            ymax = max(ymax, math.ceil(np.max(this_vals) * 1.25))
            axes[row][col].plot(this_vrest, this_vals, c=c, label=label)
            axes[row][col].set_ylabel(ordered_ylabels[i])
            axes[row][col].set_ylim([0., ymax])
            axes[row][col].set_title(ordered_subtitles[i], fontsize=mpl.rcParams['font.size'])
        if row == 1:
            axes[row][col].set_xlabel('Adjusted Vm (mV)')
            axes[row][col].set_xticks([-65, -62, -59, -56])
            axes[row][col].set_xticklabels([-65, -62, -59, -56])
    clean_axes(axes)
    fig.tight_layout()
axes[0][0].legend(loc='best', frameon=False, framealpha=0.5)
plt.show()
plt.close()



for section in raw_amplitude:
    for description in descriptions:
        plt.figure()
        for group in [group for group in raw_amplitude[section] if description in raw_amplitude[section][group]]:
            plt.scatter(rec_dict[section][group]['vrest'], raw_amplitude[section][group][description], label=group)
            plt.title('%s_%s' % (section, description))
        plt.legend(loc='best', frameon=False, framealpha=0.5)
plt.show()
plt.close()
"""

"""

ymin = {}
ymax = {}
for rec_type in rec_types:
    if rec_type in ['soma', 'dend', 'distal_spine']:
        category = 'vm'
    elif rec_type in ['soma_ina', 'dend_ina', 'ais_ina']:
        category = 'ina'
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
"""
"""
from matplotlib import cm
this_cm = cm.get_cmap()
num_items = len(rec_dict.itervalues().next()['vrest'])
colors = [this_cm(1.*i/float(num_items-1)) for i in range(num_items)]

for rec_type in rec_types:
    fig, axes = plt.subplots(2, 2, sharey=True, sharex=True)
    for i, group in enumerate([group for group in ['Control', 'TTX', 'AP5', 'AP5 + TTX'] if group in rec_dict]):
        if group not in sort_indexes:
            sort_indexes[group] = range(len(rec_dict[group]['vrest']))
            sort_indexes[group].sort(key=rec_dict[group]['vrest'].__getitem__)
            sort_indexes[group].reverse()
            rec_dict[group]['vrest'] = map(rec_dict[group]['vrest'].__getitem__, sort_indexes[group])
        if rec_type in rec_dict[group]:
            row = i / 2
            col = i % 2
            axes[row][col].set_title(group)
            rec_dict[group][rec_type] = map(rec_dict[group][rec_type].__getitem__, sort_indexes[group])
            for j, vrest in enumerate(rec_dict[group]['vrest']):
                if group == 'Control':
                    axes[row][col].plot(offset_t, rec_dict[group][rec_type][j], color=colors[j],
                                        label='%.1f mV' % vrest)
                else:
                    axes[row][col].plot(offset_t, rec_dict[group][rec_type][j], color=colors[j])
    if rec_type in ['soma', 'dend', 'distal_spine']:
        # category = 'vm'
        axes[0][0].set_ylabel('EPSP amp (mV)')
        # plt.ylim(ymin[category], ymax[category])
        plt.suptitle(rec_type+'_Vm')
    elif rec_type in ['soma_ina', 'dend_ina', 'ais_ina']:
        # category = 'ina'
        axes[0][0].set_ylabel('Current (mA/cm2)')
        # plt.ylim(ymin[category], ymax[category])
        plt.suptitle(rec_type)
    else:
        axes[0][0].set_ylabel('Current (nA)')
        plt.suptitle(rec_type)
    axes[0][0].legend(loc='best', frameon=False, framealpha=0.5)
    axes[1][0].set_xlabel('Time (ms)')
    clean_axes(axes)
    fig.tight_layout()
    fig.subplots_adjust(top=0.85, hspace=0.5, wspace=0.2)
plt.show()
plt.close()
"""
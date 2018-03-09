from specify_cells import *
import click

morph_filename = 'EB2-late-bifurcation.swc'

mech_filename = '043016 Type A - km2_NMDA_KIN5_Pr'


@click.command()
@click.option('--ttx', type=int, default=0)
@click.option('--vrest', type=float, default=-65.)
@click.option('--section', type=str, default='soma')
@click.option('--output-file-path', type=str, default=None)
def main(ttx, vrest, section, output_file_path):
    """

    :param ttx: int (bool)
    :param vrest: float
    :param section: str
    :param output_file_path: str (path)
    """

    if output_file_path is None:
        output_file_path = 'data/%s_nap_amplification_IO_DC_%s_stim_trunk.hdf5' % \
                           (datetime.datetime.today().strftime('%m%d%Y%H%M'),  section)
    globals()['output_file_path'] = output_file_path
    ttx = bool(ttx)
    globals()['ttx'] = ttx
    globals()['vrest'] = vrest
    globals()['section'] = section
    global local_random
    global i_holding
    global cell
    global rec_nodes
    global rec_locs
    global sim
    global dt
    global equilibrate
    global duration
    global v_init

    dt = 0.02
    equilibrate = 250.  # time to steady-state
    duration = 450.
    v_init = -63.
    NMDA_type = 'NMDA_KIN5'
    syn_types = ['AMPA_KIN', NMDA_type]

    local_random = random.Random()
    local_random.seed(0)

    i_holding = {'soma': 0.04, 'dend': 0.04}

    cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
    if ttx:
        # reduce_gna(0.5)
        reduce_gna(0.)
    cell.set_terminal_branch_na_gradient()

    spine_list = []
    for branch in cell.trunk:  # +cell.apical:
        spine_list.extend(branch.spines)

    # look for a trunk bifurcation
    trunk_bifurcation = [trunk for trunk in cell.trunk if cell.is_bifurcation(trunk, 'trunk')]
    if trunk_bifurcation:
        trunk_branches = [branch for branch in trunk_bifurcation[0].children if branch.type == 'trunk']
        # get where the thickest trunk branch gives rise to the tuft
        trunk = max(trunk_branches, key=lambda node: node.sec(0.).diam)
        trunk = (node for node in cell.trunk if cell.node_in_subtree(trunk, node) and 'tuft' in (child.type
                                                                                                 for child in
                                                                                                 node.children)).next()
    else:
        trunk_bifurcation = [node for node in cell.trunk if 'tuft' in (child.type for child in node.children)]
        trunk = trunk_bifurcation[0]

    spike_times = h.Vector([equilibrate + i * 10. for i in range(5)])

    stim_spines = list(spine_list)
    local_random.shuffle(stim_spines)

    rec_locs = {'soma': 0., 'dend': 0.}
    rec_nodes = {'soma': cell.tree.root, 'dend': trunk}

    for branch in cell.trunk:  # +cell.apical:
        for spine in branch.spines:
            syn = Synapse(cell, spine, syn_types, stochastic=0)
    cell.init_synaptic_mechanisms()

    sim = QuickSim(duration, verbose=True)
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)

    for description, node in rec_nodes.iteritems():
        sim.append_rec(cell, node, loc=rec_locs[description], description=description)

    sim.parameters['dt'] = dt
    sim.parameters['duration'] = duration
    sim.parameters['ttx'] = ttx
    target_vrest = offset_vm(section, vrest)
    sim.parameters['vrest'] = target_vrest

    t = np.arange(0., duration, dt)
    start = int(equilibrate/dt)
    v_th = -52.
    exceeds_th = False
    i = 0
    while not exceeds_th:
        syn = stim_spines[i].synapses[0]
        syn.source.play(spike_times)
        sim.append_rec(cell, stim_spines[i], object=stim_spines[i].synapses[0].target('AMPA_KIN'), param='_ref_g',
                       description='g_AMPA_%i' % i)
        sim.parameters['num_syns'] = i + 1
        start_time = time.time()
        sim.run(target_vrest)
        print 'Process: %i completed simulation with v_rest: %.1f, ttx: %s, num_syns: %i in %.1f s' % \
              (os.getpid(), target_vrest, str(ttx), i + 1, time.time() - start_time)
        vm = np.interp(t, sim.tvec, sim.get_rec('soma')['vec'])[start:]
        if np.max(vm) > v_th:
            exceeds_th = True
            print 'Cell crossed spike threshold with num_syns: %i' % (i + 1)
        with h5py.File(output_file_path, 'a') as f:
            sim.export_to_file(f)
        print 'EPSP amp: soma: %.2f, dend: %.2f' % (get_EPSP_amp('soma'), get_EPSP_amp('dend'))
        # sim.plot()
        i += 1


def offset_vm(description, vm_target=None):
    """

    :param description: str
    :param vm_target: float
    """
    if vm_target is None:
        vm_target = v_init
    sim.modify_stim(0, amp=0.)
    node = rec_nodes[description]
    loc = rec_locs[description]
    rec_dict = sim.get_rec(description)
    sim.modify_stim(0, node=node, loc=loc, amp=i_holding[description])
    rec = rec_dict['vec']
    offset = True
    sim.tstop = equilibrate
    t = np.arange(0., equilibrate, dt)
    sim.run(vm_target)
    vm = np.interp(t, sim.tvec, rec)
    v_rest = np.mean(vm[int((equilibrate - 3.) / dt):int((equilibrate - 1.) / dt)])
    if v_rest < vm_target - 0.5:
        i_holding[description] += 0.01
        while offset:
            if sim.verbose:
                print 'increasing i_holding to %.3f (%s)' % (i_holding[description], description)
            sim.modify_stim(0, amp=i_holding[description])
            sim.run(vm_target)
            vm = np.interp(t, sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.) / dt):int((equilibrate - 1.) / dt)])
            if v_rest < vm_target - 0.5:
                i_holding[description] += 0.01
            else:
                offset = False
    elif v_rest > vm_target + 0.5:
        i_holding[description] -= 0.01
        while offset:
            if sim.verbose:
                print 'decreasing i_holding to %.3f (%s)' % (i_holding[description], description)
            sim.modify_stim(0, amp=i_holding[description])
            sim.run(vm_target)
            vm = np.interp(t, sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.) / dt):int((equilibrate - 1.) / dt)])
            if v_rest > vm_target + 0.5:
                i_holding[description] -= 0.01
            else:
                offset = False
    sim.tstop = duration
    return v_rest


def reduce_gna(fraction):
    """

    :param fraction: float
    """
    for sec_type in ['soma', 'axon_hill', 'ais', 'axon', 'basal', 'trunk', 'apical', 'tuft']:
        for na_type in (na_type for na_type in ['nas_kin', 'nat_kin', 'nas', 'nax'] if na_type in
                cell.mech_dict[sec_type]):
            if 'value' in cell.mech_dict[sec_type][na_type]['gbar']:
                new_value = cell.mech_dict[sec_type][na_type]['gbar']['value'] * fraction
                cell.modify_mech_param(sec_type, na_type, 'gbar', new_value)
    for sec_type in ['soma', 'axon_hill', 'ais', 'axon', 'basal', 'trunk', 'apical', 'tuft']:
        for na_type in (na_type for na_type in ['nas_kin', 'nat_kin', 'nas', 'nax'] if na_type in
                cell.mech_dict[sec_type]):
            cell.reinitialize_subset_mechanisms(sec_type, na_type)


def get_EPSP_amp(description):
    """

    :param description: str
    :return: float
    """
    rec = sim.get_rec(description)['vec']
    t = np.arange(0., duration, dt)
    vm = np.interp(t, sim.tvec, rec)
    baseline = np.mean(vm[int((equilibrate - 3.) / dt):int((equilibrate - 1.) / dt)])
    peak = np.max(vm[int(equilibrate / dt):])
    return peak - baseline


if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda s: s.find(os.path.basename(__file__)) != -1, sys.argv) + 1):])
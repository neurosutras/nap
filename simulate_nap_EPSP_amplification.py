from specify_cells import *
import click

morph_filename = 'EB2-late-bifurcation.swc'

mech_filename = '043016 Type A - km2_NMDA_KIN5_Pr'


@click.command()
@click.option('--ttx', type=int, default=0)
@click.option('--ap5', type=int, default=0)
@click.option('--section', type=str, default='soma')
@click.option('--output-file-path', type=str, default=None)
def main(ttx, ap5, section, output_file_path):
    """

    :param ttx: int (bool)
    :param ap5: int (bool)
    :param section: str
    :param output_file_path: str (path)
    """

    if output_file_path is None:
        output_file_path = 'data/%s_nap_amplification_DC_%s_stim_trunk.hdf5' % \
                           (datetime.datetime.today().strftime('%m%d%Y%H%M'), section)
    globals()['output_file_path'] = output_file_path
    ttx = bool(ttx)
    ap5 = bool(ap5)
    globals()['ttx'] = ttx
    globals()['ap5'] = ap5
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

    if not ap5:
        syn_types = ['AMPA_KIN', NMDA_type]
    else:
        syn_types = ['AMPA_KIN']

    local_random = random.Random()
    local_random.seed(0)

    i_holding = {'soma': 0.04, 'dend': 0.04}

    cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
    if ttx:
        reduce_gna(0.5)
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

    num_syns_dict = {'soma': {True: 7, False: 5}, 'dend': {True: 4, False: 3}}
    num_syns = num_syns_dict[section][ap5]
    # min_num_syns = np.min([val for element in num_syns_dict.itervalues() for val in element.itervalues()])

    stim_spines = local_random.sample(spine_list, num_syns)
    """
    stim_spine_locs = np.array([abs(cell.get_distance_to_node(cell.tree.root, trunk) -
                                    cell.get_distance_to_node(cell.tree.root, spine)) for spine in stim_spines])
    index = np.where(stim_spine_locs[:min_num_syns] == np.max(stim_spine_locs[:min_num_syns]))[0][0]
    distal_spine = stim_spines[index]
    """

    rec_locs = {'soma': 0., 'dend': 0., 'spine': 0.5}
    rec_nodes = {'soma': cell.tree.root, 'dend': trunk}

    for branch in cell.trunk:  # +cell.apical:
        for spine in branch.spines:
            syn = Synapse(cell, spine, syn_types, stochastic=0)
    cell.init_synaptic_mechanisms()

    sim = QuickSim(duration, verbose=True)
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)

    for description, node in rec_nodes.iteritems():
        sim.append_rec(cell, node, loc=rec_locs[description], description=description)

    sim.append_rec(cell, cell.tree.root, object=cell.tree.root.sec(0.5), param='_ref_ina_nas', description='soma_ina')
    sim.append_rec(cell, cell.axon[1], object=cell.axon[1].sec(0.5), param='_ref_ina_nax', description='ais_ina')
    sim.append_rec(cell, trunk, object=trunk.sec(0.5), param='_ref_ina_nas', description='dend_ina')

    for i, spine in enumerate(stim_spines):
        syn = spine.synapses[0]
        syn.source.play(spike_times)
        sim.append_rec(cell, spine, loc=rec_locs['spine'], description='spine_%i' % i)
        if not ap5:
            sim.append_rec(cell, spine, object=syn.target(NMDA_type), param='_ref_i', description='i_NMDA_%i' % i)

    sim.parameters['dt'] = dt
    sim.parameters['duration'] = duration
    sim.parameters['ap5'] = ap5
    sim.parameters['ttx'] = ttx
    sim.parameters['DC_section'] = section

    """
    delta_vm = 0.  # 8.  # 0.
    start_time = time.time()
    target_vrest = v_init + delta_vm
    target_vrest = offset_vm(section, target_vrest)
    sim.run(target_vrest)

    print 'Process: %i completed simulation with v_rest: %.1f, ap5: %s, ttx: %s in %.i s' % (os.getpid(), target_vrest,
                                                                                             str(ap5), str(ttx),
                                                                                             time.time() - start_time)
    print 'EPSP amp: soma: %.2f, dend: %.2f' % (get_EPSP_amp('soma'), get_EPSP_amp('dend'))
    sim.plot()
    """

    for delta_vm in np.arange(-3., 8., 2.):
        start_time = time.time()
        target_vrest = v_init + delta_vm
        target_vrest = offset_vm(section, target_vrest)
        sim.run(target_vrest)
        sim.parameters['vrest'] = target_vrest
        with h5py.File(output_file_path, 'a') as f:
            sim.export_to_file(f)
        print 'Process: %i completed simulation with v_rest: %.1f, ap5: %s, ttx: %s in %.i s' % \
              (os.getpid(), target_vrest, str(ap5), str(ttx), time.time() - start_time)
        print 'EPSP amp: soma: %.2f, dend: %.2f' % (get_EPSP_amp('soma'), get_EPSP_amp('dend'))
        # sim.plot()


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
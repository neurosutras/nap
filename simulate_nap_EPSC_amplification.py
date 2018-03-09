from specify_cells import *
import click

morph_filename = 'EB2-late-bifurcation.swc'

mech_filename = '043016 Type A - km2_NMDA_KIN5_Pr'

context = Context()


@click.command()
@click.option('--ttx', type=int, default=0)
@click.option('--bar', type=int, default=0)
@click.option('--zd', type=int, default=0)
@click.option('--section', type=str, default='soma')
@click.option('--output-file-path', type=str, default=None)
@click.option('--tune-imax', is_flag=True)
@click.option('--debug', is_flag=True)
@click.option('--plot', is_flag=True)
def main(ttx, bar, zd, section, output_file_path, tune_imax, debug, plot):
    """

    :param ttx: int (bool)
    :param bar: int (bool)
    :param zd: int (bool)
    :param section: str
    :param output_file_path: str (path)
    :param tune_imax: bool
    :param debug: bool
    :param plot: bool
    """

    if output_file_path is None:
        output_file_path = 'data/%s_nap_EPSC_amplification_DC_%s_stim_soma.hdf5' % \
                           (datetime.datetime.today().strftime('%m%d%Y%H%M'), section)
    ttx = bool(ttx)
    bar = bool(bar)
    zd = bool(zd)

    dt = 0.02
    equilibrate = 250.  # time to steady-state
    duration = 450.
    v_init = -63.
    syn_types = ['EPSC']

    local_random = random.Random()
    local_random.seed(0)

    i_holding = {'soma': 0.04, 'dend': 0.04}

    cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

    cell.set_terminal_branch_na_gradient()
    if ttx:
        multiply_mech_by_fraction(cell, ['nas_kin', 'nat_kin', 'nas', 'nax'], 'gbar', 0.)
    if bar:
        multiply_mech_by_fraction(cell, ['kap', 'kad'], 'gkabar', 0.65)
    if zd:
        multiply_mech_by_fraction(cell, ['h'], 'ghbar', 0.6)

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

    rec_locs = {'soma': 0.}
    rec_nodes = {'soma': cell.tree.root}

    syn = Synapse(cell, cell.tree.root, syn_types, stochastic=False)

    sim = QuickSim(duration, verbose=True)
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)

    for description, node in rec_nodes.iteritems():
        sim.append_rec(cell, node, loc=rec_locs[description], description=description)

    # sim.append_rec(cell, cell.tree.root, object=cell.tree.root.sec(0.5), param='_ref_ina_nas', description='soma_ina')
    sim.append_rec(cell, cell.axon[1], object=cell.axon[1].sec(0.5), param='_ref_ina_nax', description='ais_ina')
    sim.append_rec(cell, cell.tree.root, object=cell.tree.root.sec(0.5), param='_ref_ik_kap', description='soma_ika')
    sim.append_rec(cell, cell.tree.root, object=cell.tree.root.sec(0.5), param='_ref_i_h', description='soma_ih')

    syn.source.play(spike_times)

    sim.parameters['dt'] = dt
    sim.parameters['duration'] = duration
    sim.parameters['ttx'] = ttx
    sim.parameters['bar'] = bar
    sim.parameters['zd'] = zd
    sim.parameters['DC_section'] = section

    context.update(locals())

    if tune_imax:
        start_time = time.time()
        target_vrest = offset_vm(section, v_init)
        x0 = [0.05]
        result = optimize.minimize(EPSP_amp_error, x0, method='Nelder-Mead', options={'fatol': 1e-3, 'xatol': 1e-3,
                                                                                      'disp': True, 'maxiter': 20})
        print result
        print 'Process: %i completed imax optimization with imax: %.4f, v_rest: %.1f, ttx: %s, bar: %s, zd: %s in %.i' \
              ' s' % (os.getpid(), result.x[0], target_vrest, str(ttx), str(bar), str(zd), time.time() - start_time)
    else:
        if bar:
            if ttx:
                imax = 0.0493  # 0.059765625
            else:
                imax = 0.0462  # 0.053750
        elif zd:
            if ttx:
                imax = 0.0545  # 0.0803125
            else:
                imax = 0.0513  # 0.076250
        elif ttx:
            imax = 0.0555  # 0.08609375
        else:
            imax = 0.0523  # 0.081875
        syn.target('EPSC').imax = imax

        if debug:
            delta_vm = 0.
            start_time = time.time()
            target_vrest = v_init + delta_vm
            target_vrest = offset_vm(section, target_vrest)
            sim.run(target_vrest)
            print 'Process: %i completed simulation with v_rest: %.1f, ttx: %s, bar: %s, zd: %s in %.i s' % \
                  (os.getpid(), target_vrest, str(ttx), str(bar), str(zd), time.time() - start_time)
            print 'EPSP amp: soma: %.2f' % (get_EPSP_amp('soma'))
            if plot:
                sim.plot()
        else:
            for delta_vm in np.arange(-3., 8., 2.): # [6.]:
                start_time = time.time()
                target_vrest = v_init + delta_vm
                target_vrest = offset_vm(section, target_vrest)
                sim.run(target_vrest)
                sim.parameters['vrest'] = target_vrest
                with h5py.File(output_file_path, 'a') as f:
                    sim.export_to_file(f)
                print 'Process: %i completed simulation with v_rest: %.1f, ttx: %s, bar: %s, zd: %s in %.i s' % \
                  (os.getpid(), target_vrest, str(ttx), str(bar), str(zd), time.time() - start_time)
                print 'EPSP amp: soma: %.3f' % (get_EPSP_amp('soma'))
                if plot:
                    sim.plot()

    context.update(locals())


def offset_vm(description, vm_target=None):
    """

    :param description: str
    :param vm_target: float
    """
    if vm_target is None:
        vm_target = context.v_init
        context.sim.modify_stim(0, amp=0.)
    node = context.rec_nodes[description]
    loc = context.rec_locs[description]
    rec_dict = context.sim.get_rec(description)
    context.sim.modify_stim(0, node=node, loc=loc, amp=context.i_holding[description])
    rec = rec_dict['vec']
    offset = True
    context.sim.tstop = context.equilibrate
    t = np.arange(0., context.equilibrate, context.dt)
    context.sim.run(vm_target)
    vm = np.interp(t, context.sim.tvec, rec)
    v_rest = np.mean(vm[int((context.equilibrate - 3.) / context.dt):int((context.equilibrate - 1.) / context.dt)])
    if v_rest < vm_target - 0.5:
        context.i_holding[description] += 0.01
        while offset:
            if context.sim.verbose:
                print 'increasing i_holding to %.3f (%s)' % (context.i_holding[description], description)
            context.sim.modify_stim(0, amp=context.i_holding[description])
            context.sim.run(vm_target)
            vm = np.interp(t, context.sim.tvec, rec)
            v_rest = np.mean(vm[int((context.equilibrate - 3.) / context.dt):int((context.equilibrate - 1.) /
                                                                                 context.dt)])
            if v_rest < vm_target - 0.5:
                context.i_holding[description] += 0.01
            else:
                offset = False
    elif v_rest > vm_target + 0.5:
        context.i_holding[description] -= 0.01
        while offset:
            if context.sim.verbose:
                print 'decreasing i_holding to %.3f (%s)' % (context.i_holding[description], description)
            context.sim.modify_stim(0, amp=context.i_holding[description])
            context.sim.run(vm_target)
            vm = np.interp(t, context.sim.tvec, rec)
            v_rest = np.mean(vm[int((context.equilibrate - 3.) / context.dt):int((context.equilibrate - 1.) /
                                                                                 context.dt)])
            if v_rest > vm_target + 0.5:
                context.i_holding[description] -= 0.01
            else:
                offset = False
    context.sim.tstop = context.duration
    return v_rest


def multiply_mech_by_fraction(cell, mech_name_list, param_name, fraction):
    """

    :param cell:
    :param mech_name_list:
    :param param_name:
    :param fraction:
    """
    for sec_type in ['soma', 'axon_hill', 'ais', 'axon', 'basal', 'trunk', 'apical', 'tuft']:
        for mech_name in (mech_name for mech_name in mech_name_list if mech_name in cell.mech_dict[sec_type]):
            for node in cell.get_nodes_of_subtype(sec_type):
                for seg in node.sec:
                    if hasattr(seg, mech_name):
                        mech = getattr(seg, mech_name)
                        if hasattr(mech, param_name):
                            val = getattr(mech, param_name)
                            setattr(mech, param_name, val * fraction)


def get_EPSP_amp(description):
    """

    :param description: str
    :return: float
    """
    rec = context.sim.get_rec(description)['vec']
    t = np.arange(0., context.duration, context.dt)
    vm = np.interp(t, context.sim.tvec, rec)
    baseline = np.mean(vm[int((context.equilibrate - 3.) / context.dt):int((context.equilibrate - 1.) / context.dt)])
    peak = np.max(vm[int(context.equilibrate / context.dt):])
    return peak - baseline


def EPSP_amp_error(x):
    """

    :param x: array
    :return: float
    """
    start_time = time.time()
    context.syn.target('EPSC').imax = x[0]
    context.sim.run(context.v_init)
    EPSP_amp = get_EPSP_amp('soma')
    Err = ((EPSP_amp - 2.) / 0.1) ** 2.
    print 'imax: %.4f; EPSP_amp: %.3f; error: %.4E' % (x[0], EPSP_amp, Err)
    return Err


if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda s: s.find(os.path.basename(__file__)) != -1, sys.argv) + 1):])
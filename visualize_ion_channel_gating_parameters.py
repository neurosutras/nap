__author__ = 'milsteina'
from plot_results import *


def visualize_kap(new_params=None):
    """
    
    :param new_params:  
    """
    current = {'celsius': 35., 'temp': 24., 'vhalfn': 11., 'a0n': 0.05, 'zetan': -1.5, 'gmn': 0.55, 'pw': -1.,
               'tq': -40., 'qq': 5., 'nmin': 0.1, 'nscale': 1., 'vhalfl': -56., 'a0l': 0.05, 'zetal': 3., 'lmin': 2.,
               'lscale': 1.,'q10': 5.}
    if new_params is not None:
        for param in current:
            if param in new_params:
                current[param] = new_params[param]
    v = np.arange(-100., 60., 1.)
    qt = current['q10'] ** ((current['celsius'] - current['temp']) / 10.)
    zeta = current['zetan'] + current['pw'] / (1. + np.exp((v - current['tq']) / current['qq']))
    alpn = np.exp(zeta * (v - current['vhalfn']) * 1.e-3 * 9.648e4 / (8.315 * (273.16 + current['celsius'])))
    betn = np.exp(zeta * current['gmn'] * (v - current['vhalfn']) * 1.e-3 * 9.648e4 /
                  (8.315 * (273.16 + current['celsius'])))
    an = alpn
    ninf = 1. / (1. + an)
    taun = betn / (qt * current['a0n'] * (1. + an))
    taun_min_indexes = np.where(taun < current['nmin'])[0]
    taun[taun_min_indexes] = current['nmin']
    taun /= current['nscale']

    alpl = np.exp(current['zetal'] * (v - current['vhalfl']) * 1.e-3 * 9.648e4 /
                  (8.315 * (273.16 + current['celsius'])))
    al = alpl
    linf = 1. / (1. + al)
    taul = 0.26 * (v + 50.)
    taul_min_indexes = np.where(taul < current['lmin'])[0]
    taul[taul_min_indexes] = current['lmin']
    taul /= current['lscale']
    fig, axes = plt.subplots(1, 2)
    axes[0].plot(v, ninf, label='Activation gate')
    axes[0].plot(v, linf, label='Inactivation gate')
    axes[1].plot(v, taun, label='Activation tau')
    axes[1].plot(v, taul, label='Inactivation tau')
    axes[0].set_ylabel('Normalized conductance')
    axes[1].set_ylabel('Time constant (ms)')
    axes[0].set_xlabel('Voltage (mV)')
    axes[1].set_xlabel('Voltage (mV)')
    fig.suptitle('Ka')
    axes[0].legend(loc='best', frameon=False, framealpha=0.5)
    axes[1].legend(loc='best', frameon=False, framealpha=0.5)
    clean_axes(axes)
    fig.tight_layout()
    plt.show()
    plt.close()


def visualize_h(new_params=None):
    """

    :param new_params:
    """
    current = {'celsius': 35., 'eh': 30., 'vhalfl': -90., 'vhalft': -75., 'a0t': 0.011, 'zetal': 4., 'zetat': 2.2,
               'gmt': 0.4, 'q10': 4.5, 'qtl': 1.}
    if new_params is not None:
        for param in current:
            if param in new_params:
                current[param] = new_params[param]
    v = np.arange(-120., 60., 1.)
    qt = current['q10'] ** ((current['celsius'] - 33.) / 10.)
    alpt = np.exp(1e-3 * current['zetat'] * (v - current['vhalft']) * 9.648e4 / (8.315 * (273.16 + current['celsius'])))
    alpl = np.exp(1e-3 * current['zetal'] * (v - current['vhalfl']) * 9.648e4 / (8.315 * (273.16 + current['celsius'])))
    bett = np.exp(1e-3 * current['zetat'] * current['gmt'] * (v - current['vhalft']) * 9.648e4 /
                  (8.315 * (273.16 + current['celsius'])))
    a = alpt
    linf = 1. / (1. + alpl)
    taul = bett / (current['qtl'] * qt * current['a0t'] * (1. + a))

    fig, axes = plt.subplots(1, 2)
    axes[0].plot(v, linf, label='Activation gate')
    axes[1].plot(v, taul, label='Aactivation tau')
    axes[0].set_ylabel('Normalized conductance')
    axes[1].set_ylabel('Time constant (ms)')
    axes[0].set_xlabel('Voltage (mV)')
    axes[1].set_xlabel('Voltage (mV)')
    fig.suptitle('h')
    axes[0].legend(loc='best', frameon=False, framealpha=0.5)
    axes[1].legend(loc='best', frameon=False, framealpha=0.5)
    clean_axes(axes)
    fig.tight_layout()
    plt.show()
    plt.close()


def visualize_kdr(new_params=None):
    """
    
    :param new_params:  
    """
    orig = {'celsius': 35., 'temp': 24., 'vhalfn': 13., 'a0n': 0.02, 'zetan': -3., 'gmn': 0.7, 'nmin': 2.}
    current = {'celsius': 35., 'temp': 24., 'vhalfn': -10., 'a0n': 0.055, 'zetan': -5., 'gmn': 0.9, 'nmin': 5.}
    if new_params is not None:
        for param in current:
            if param in new_params:
                current[param] = new_params[param]
    v = np.arange(-80., 40., 1.)
    alpn = lambda v, zetan, vhalfn, celsius: np.exp(zetan * (v - vhalfn) * 1.e-3 * 9.648e4 / (8.315 * (273.16 +
                                                                                                       celsius)))
    betn = lambda v, zetan, gmn, vhalfn, celsius: np.exp(zetan * gmn * (v - vhalfn) * 1.e-3 * 9.648e4 / (8.315 *
                                                                                                (273.16 + celsius)))
    fig, axes = plt.subplots(1, 2)
    for param_set, set_label in zip([orig, current], ['orig', 'new']):
        a = alpn(v, param_set['zetan'], param_set['vhalfn'], param_set['celsius'])
        ninf = 1. / (1. + a)
        taun = param_set['nmin'] + betn(v, param_set['zetan'], param_set['gmn'], param_set['vhalfn'], param_set['celsius']) / \
               (param_set['a0n'] * (1. + a))
        axes[0].plot(v, ninf, label=set_label)
        axes[1].plot(v, taun, label=set_label)
    axes[0].set_title('ninf')
    axes[1].set_title('ntau')
    plt.legend(loc='best', frameon=False, framealpha=0.5)
    plt.show()
    plt.close()


def visualize_nax(new_params=None):
    """

    """
    orig = {'celsius': 35., 'tha': -30., 'qa': 7.2, 'Ra': 0.4, 'Rb': 0.124, 'thi1': -45., 'thi2': -45., 'qd': 1.5,
            'qg': 1.5, 'mmin': 0.02, 'hmin': 0.5, 'q10': 2., 'Rg': 0.01, 'Rd': 0.03, 'thinf': -50., 'qinf': 4.,
            'sh': 0., 'sha': 0.}
    current = {'celsius': 35., 'tha': -30., 'qa': 7.2, 'Ra': 0.4, 'Rb': 0.124, 'thi1': -45., 'thi2': -45., 'qd': 1.5,
               'qg': 1.5, 'mmin': 0.02, 'hmin': 0.5, 'q10': 2., 'Rg': 0.01, 'Rd': 0.03, 'thinf': -50., 'qinf': 4.,
               'sh': 5.514521615045231, 'sha': 0.}
    if new_params is not None:
        for param in current:
            if param in new_params:
                current[param] = new_params[param]
    v = np.arange(-80., 40., 1.)
    trap0 = lambda v, th, R, q: R * (v - th) / (1. - np.exp(-(v - th) / q)) if abs(v - th) > 1.e-6 else R * q
    trap0 = np.vectorize(trap0, excluded=[1,2,3])

    fig, axes = plt.subplots(2, 2)
    for param_set, set_label in zip([orig, current], ['orig', 'new']):
        qt = param_set['q10'] ** ((param_set['celsius'] - 24.) / 10.)
        a = trap0(v, param_set['tha'] + param_set['sh'] + param_set['sha'], param_set['Ra'], param_set['qa'])
        b = trap0(-v, -(param_set['tha'] + param_set['sh'] + param_set['sha']), param_set['Rb'], param_set['qa'])
        mtau = 1. / (a + b) / qt
        indexes = np.where(mtau < param_set['mmin'])[0]
        if np.any(indexes):
            mtau[indexes] = param_set['mmin']
        minf = a / (a + b)
        axes[0][0].plot(v, minf, label=set_label+': minf')
        axes[1][0].plot(v, mtau, label=set_label+': mtau')
        a = trap0(v, param_set['thi1'] + param_set['sh'], param_set['Rd'], param_set['qd'])
        b = trap0(-v, -(param_set['thi2'] + param_set['sh']), param_set['Rg'], param_set['qg'])
        htau = 1. / (a + b) / qt
        indexes = np.where(htau < param_set['hmin'])[0]
        if np.any(indexes):
            htau[indexes] = param_set['hmin']
        hinf = 1. / (1. + np.exp((v - param_set['thinf'] - param_set['sh']) / param_set['qinf']))
        axes[0][0].plot(v, hinf, label=set_label+': hinf', linestyle='--')
        axes[1][1].plot(v, htau, label=set_label+': htau', linestyle='--')
    axes[0][0].legend(loc='best', frameon=False, framealpha=0.5)
    axes[1][0].legend(loc='best', frameon=False, framealpha=0.5)
    axes[1][1].legend(loc='best', frameon=False, framealpha=0.5)
    clean_axes(axes)
    fig.tight_layout()
    plt.show()
    plt.close()


def visualize_km2(new_params=None):
    """

    """
    orig = {'vhalfl': -42., 'kl': -4.,'vhalft': -42., 'a0t': 0.009, 'zetat': 7., 'gmt': .4,
            't0': 15., 'st': 1.}
    current = {'vhalfl': -40., 'kl': -7., 'vhalft': -42., 'a0t': 0.009, 'zetat': 7., 'gmt': .4,
               't0': 15., 'st': 1.}
    if new_params is not None:
        for param in current:
            if param in new_params:
                current[param] = new_params[param]
    v = np.arange(-80., 40., 1.)
    alpt = lambda v, vhalft, zetat: np.exp(0.0378 * zetat * (v - vhalft))
    bett = lambda v, vhalft, zetat, gmt: np.exp(0.0378 * zetat * gmt * (v - vhalft))
    fig, axes = plt.subplots(2, 1)
    for param_set, set_label in zip([orig, current], ['orig', 'new']):
        inf = (1. / (1. + np.exp((v - param_set['vhalfl']) / param_set['kl'])))
        a = alpt(v, param_set['vhalft'], param_set['zetat'])
        tau = param_set['t0'] + bett(v, param_set['vhalft'], param_set['zetat'], param_set['gmt']) / \
                                (param_set['a0t'] * (1. + a))
        axes[0].plot(v, inf, label=set_label)
        axes[1].plot(v, tau, label=set_label)
    axes[0].set_title('inf')
    axes[1].set_title('tau')
    plt.legend(loc='best', frameon=False, framealpha=0.5)
    plt.show()
    plt.close()


def visualize_CadepK():
    B = .26
    gbkbar = .0003
    gskbar = .0005
    gcakmult = 1.
    ca0 = .00007
    tau = 9.
    taucadiv = 1.
    tauskdiv = 1.
    ask = 1.
    bsk = 1.
    alphar = 7.5
    stau = 10.
    v = np.arange(-80., 40., 1.)
    ca_i = np.arange(0.00005, 0.1, 0.00001)

    exp1 = lambda A, d, k, x: A/np.exp((12.*np.log10(x)+d)/k)
    alphaq = lambda x: exp1(0.00246, 28.48, -4.5, x)
    betaq = lambda x: exp1(0.006, 60.4, 35, x)
    plt.figure()
    plt.plot(ca_i, alphaq(ca_i))
    plt.figure()
    plt.plot(ca_i, betaq(ca_i))

    betar = lambda v: 0.11 / np.exp((v - 35.) / 14.9)
    plt.figure()
    plt.plot(v, betar(v))

    sinf = lambda x: 1. / (1. + 4. / (1000. * x))
    plt.figure()
    plt.plot(ca_i, sinf(ca_i))

    plt.show()
    plt.close()


def visualize_NMDAR_g_V(Kd=9.98, gamma=0.101, mg=1.):
    """

    :param Kd: float
    :param gamma: float
    :param mg:  float
    """
    v = np.arange(-100., 50., 1.)
    B = 1. / (1. + np.exp(gamma * (-v)) * (mg / Kd))
    B /= np.max(B)
    fig, axes = plt.subplots(1)
    axes.plot(v, B)
    axes.set_ylabel('Normalized conductance')
    axes.set_xlabel('Voltage (mV)')
    axes.set_title('NMDAR g-V')
    clean_axes(axes)
    plt.show()
    plt.close()


def visualize_nax_gates():
    """

    """
    param_set = {'celsius': 35., 'tha': -30., 'qa': 7.2, 'Ra': 0.4, 'Rb': 0.124, 'thi1': -45., 'thi2': -45., 'qd': 1.5,
                 'qg': 1.5, 'mmin': 0.02, 'hmin': 0.5, 'q10': 2., 'Rg': 0.01, 'Rd': 0.03, 'thinf': -50., 'qinf': 4.,
                 'sh': 1.7}
    v = np.arange(-100., 50., 1.)
    trap0 = np.vectorize(
        lambda v, th, R, q: R * (v - th) / (1. - np.exp(-(v - th) / q)) if abs(v - th) > 1.e-6 else R * q,
        excluded=range(1,4))

    fig, axes = plt.subplots(1)
    qt = param_set['q10'] ** ((param_set['celsius'] - 24.) / 10.)
    a = trap0(v, param_set['tha']+param_set['sh'], param_set['Ra'], param_set['qa'])
    b = trap0(-v, -(param_set['tha']+param_set['sh']), param_set['Rb'], param_set['qa'])
    mtau = 1. / (a + b) / qt
    minf = a / (a + b)
    axes.plot(v, minf, label='m (activation)')
    a = trap0(v, param_set['thi1']+param_set['sh'], param_set['Rd'], param_set['qd'])
    b = trap0(-v, -(param_set['thi2']+param_set['sh']), param_set['Rg'], param_set['qg'])
    htau = 1. / (a + b) / qt
    hinf = 1. / (1. + np.exp((v - (param_set['thinf']+param_set['sh'])) / param_set['qinf']))
    axes.plot(v, hinf, label='h (inactivation)')
    axes.set_ylabel('Normalized conductance')
    axes.set_xlabel('Voltage (mV)')
    axes.set_title('Na channel g-V')
    clean_axes(axes)
    plt.legend(loc='best', frameon=False, framealpha=0.5)
    plt.show()
    plt.close()
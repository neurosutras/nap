__author__ = 'milsteina'
try:
    from mpi4py import MPI
except Exception:
    pass
import h5py
import math
import pickle
import datetime
import copy
import time
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
import scipy.signal as signal
import random
import pprint
import sys
import os


#---------------------------------------Some global variables and functions------------------------------

data_dir = 'data/'
morph_dir = 'morphologies/'

freq = 100      # Hz, frequency at which AC length constant will be computed
d_lambda = 0.1  # no segment will be longer than this fraction of the AC length constant

"""
Structure of Mechanism Dictionary: dict of dicts

keys:               description:
'mechanism name':   Value is dictionary specifying how to set parameters at the mechanism level.
'cable':            Value is dictionary specifying how to set basic cable parameters at the section level. Includes
                        'Ra', 'cm', and the special parameter 'spatial_res', which scales the number of segments per
                        section for the specified sec_type by a factor of an exponent of 3.
'ions':             Value is dictionary specifying how to set parameters for ions at the section or segment level.
                    These parameters must be specified **after** all other mechanisms have been inserted.
values:
None:               Use default values for all parameters within this mechanism.
dict:
    keys:
    'parameter name':
    values:     dict:
                        keys:        value:
                        'origin':   'self':     Use 'value' as a baseline value.
                                    sec_type:   Inherit value from last seg of the closest node with sec of
                                                sec_type along the path to root.
                        'value':    float:      If 'origin' is 'self', contains the baseline value.
                        'slope':    float:      If exists, contains slope in units per um. If not, use
                                                constant 'value' for the full length of sec.
                        'max':      float:      If 'slope' exists, 'max' is an upper limit for the value
                        'min':      float:      If 'slope' exists, min is a lower limit for the value

"""

default_mech_dict = {'ais': {'cable': {'Ra': {'origin': 'soma'}, 'cm': {'origin': 'soma'}},
                             'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}},
                     'apical': {'cable': {'Ra': {'origin': 'soma'}, 'cm': {'origin': 'soma'}},
                                'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}},
                     'axon': {'cable': {'Ra': {'origin': 'soma'}, 'cm': {'origin': 'soma'}},
                              'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}},
                     'axon_hill': {'cable': {'Ra': {'origin': 'soma'}, 'cm': {'origin': 'soma'}},
                              'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}},
                     'basal': {'cable': {'Ra': {'origin': 'soma'}, 'cm': {'origin': 'soma'}},
                               'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}},
                     'soma': {'cable': {'Ra': {'value': 150.}, 'cm': {'value': 1.}},
                              'pas': {'e': {'value': -67.}, 'g': {'value': 2.5e-05}}},
                     'trunk': {'cable': {'Ra': {'origin': 'soma'}, 'cm': {'origin': 'soma'}},
                               'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}},
                     'tuft': {'cable': {'Ra': {'origin': 'soma'}, 'cm': {'origin': 'soma'}},
                              'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}},
                     'spine_neck': {'cable': {'Ra': {'origin': 'soma'}, 'cm': {'origin': 'soma'}},
                              'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}},
                     'spine_head': {'cable': {'Ra': {'origin': 'soma'}, 'cm': {'origin': 'soma'}},
                              'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}}}


def lambda_f(sec, f=freq):
    """
    Calculates the AC length constant for the given section at the frequency f
    Used to determine the number of segments per hoc section to achieve the desired spatial and temporal resolution
    :param sec : :class:'h.Section'
    :param f : int
    :return : int
    """
    diam = np.mean([seg.diam for seg in sec])
    Ra = sec.Ra
    cm = np.mean([seg.cm for seg in sec])
    return 1e5*math.sqrt(diam/(4.*math.pi*f*Ra*cm))


def d_lambda_nseg(sec, lam=d_lambda, f=freq):
    """
    The AC length constant for this section and the user-defined fraction is used to determine the maximum size of each
    segment to achieve the d esired spatial and temporal resolution. This method returns the number of segments to set
    the nseg parameter for this section. For tapered cylindrical sections, the diam parameter will need to be
    reinitialized after nseg changes.
    :param sec : :class:'h.Section'
    :param lam : int
    :param f : int
    :return : int
    """
    L = sec.L
    return int((L/(lam*lambda_f(sec, f))+0.9)/2)*2+1


def scaleSWC(filenameBase, mag=100, scope='neurolucida'):
    # this function rescales the SWC file with the real distances.
    f = open(morph_dir+filenameBase+'.swc')
    lines = f.readlines()
    f.close()
    Points = []
    if mag == 100:
        if scope == 'neurolucida':
            xyDist = 0.036909375  # 0.07381875
            zDist = 1.0
        else:
            xyDist = 0.065
            zDist = 0.05
    else:
        raise Exception('Calibration for {}X objective unknown.'.format(mag))
    for line in lines:
        ll = line.split(' ')
        nn = int(float(ll[0]))    # label of the point
        tp = int(float(ll[1]))  # point type
        py = float(ll[2])    # note the inversion of x, y.
        px = float(ll[3])
        z = float(ll[4])    # z
        r = float(ll[5])    # radius of the sphere.
        np = int(float(ll[6]))    # parent point id.
        # get the length in micron
        py *= xyDist; px *= xyDist; r = r*xyDist; z *= zDist
        Points.append([nn,tp,py,px,z,r,np])

    print 'Saving SWC to file '+filenameBase+'-scaled.swc'
    f = open(morph_dir+filenameBase+'-scaled.swc', 'w')
    for [nn,tp,py,px,z,r,np] in Points:
        ll = str(int(nn))+' '+str(int(tp))+' '+str(py)+' '+str(px)+' '+str(z)+' '+str(r)+' '+str(int(np))+'\n'
        f.write(ll)
    f.close()


def investigateSWC(filenameBase):
    # this function reports the min and max values for y, x, z, and radius from an SWC file.
    f = open(morph_dir+filenameBase+'.swc')
    lines = f.readlines()
    f.close()
    xvals = []
    yvals = []
    zvals = []
    rvals = []
    for line in lines:
        ll = line.split(' ')
        yvals.append(float(ll[2]))    # note the inversion of x, y.
        xvals.append(float(ll[3]))
        zvals.append(float(ll[4]))    # z
        rvals.append(float(ll[5]))    # radius of the sphere.
    print 'x - ',min(xvals),':',max(xvals)
    print 'y - ',min(yvals),':',max(yvals)
    print 'z - ',min(zvals),':',max(zvals)
    print 'r - ',min(rvals),':',max(rvals)


def translateSWCs():
    """
    Eric Bloss has produced high resolution .swc files that each contain a volume 10 um deep in z. This method
    determines from the filename the z offset of each file and translates the z coordinates of the .swc files to
    facilitate stitching them together into a single volume. Also changes the sec_type of any node that is not a root
    and has no children within a file to 7 to indicate a leaf that potentially needs to be connected to a nearby root.
    Also attempts to connect unconnected nodes that are within 2 um of each other across consecutive slices, and labels
    them with sec_type = 8. This doesn't work particularly well and files must be extensively proofread in NeuTuMac.
    """
    num_nodes = 0
    outputname = 'combined-offset-connected.swc'
    out_f = open(outputname, 'w')
    # out_test = open('combined-offset-connected.swc', 'w')
    prev_nodes = {}
    filenames = []
    z_offsets = []
    for filename in os.listdir('.'):
        if '.swc' in filename and not '-offset' in filename:
            filenames.append(filename)
            z_offsets.append(float(filename.split('z=')[1].split(' ')[0])/10.0)
    indexes = range(len(z_offsets))
    indexes.sort(key=z_offsets.__getitem__)
    for i in indexes:
        f = open(filenames[i])
        lines = f.readlines()
        f.close()
        num_nodes += len(prev_nodes)
        nodes = {}
        leaves = []
        for line in [line.split(' ') for line in lines if not line.split(' ')[0] in ['#', '\r\n']]:
            index = int(float(line[0])) + num_nodes  # node index
            nodes[index] = {}
            nodes[index]['type'] = int(float(line[1]))  # sec_type
            nodes[index]['y'] = float(line[2])  # note the inversion of x, y.
            nodes[index]['x'] = float(line[3])
            nodes[index]['z'] = float(line[4]) + z_offsets[i]
            nodes[index]['r'] = float(line[5])  # radius of the sphere.
            nodes[index]['parent'] = int(float(line[6]))  # index of parent node
            if not nodes[index]['parent'] == -1:
                nodes[index]['parent'] += num_nodes
                leaves.append(index)
        for index in nodes:  # keep nodes with no children
            parent = nodes[index]['parent']
            if parent in leaves:
                leaves.remove(parent)
        for index in leaves:
            nodes[index]['type'] = 7
        print 'Saving '+filenames[i]+' to '+outputname
        if prev_nodes:
            leaves = [index for index in nodes if (nodes[index]['type'] == 7 or nodes[index]['parent'] == -1)]
            for prev_index in [index for index in prev_nodes if (prev_nodes[index]['type'] == 7 or
                                                                prev_nodes[index]['parent'] == -1)]:
                for index in leaves:
                    distance = math.sqrt((prev_nodes[prev_index]['x']-nodes[index]['x'])**2 +
                                         (prev_nodes[prev_index]['y']-nodes[index]['y'])**2 +
                                         (prev_nodes[prev_index]['z']-nodes[index]['z'])**2)
                    # print prev_index, index, distance
                    if distance < 2.:
                        prev_nodes[prev_index]['type'] = 8
                        nodes[index]['type'] = 8
                        nodes[index]['parent'] = prev_index
                        leaves.remove(index)
                        break
        for index in prev_nodes:
            line = str(index)+' '+str(prev_nodes[index]['type'])+' '+str(prev_nodes[index]['y'])+' '+\
                   str(prev_nodes[index]['x'])+' '+str(prev_nodes[index]['z'])+' '+str(prev_nodes[index]['r'])+' '+\
                   str(prev_nodes[index]['parent'])+'\n'
            out_f.write(line)
        prev_nodes = copy.deepcopy(nodes)
    for index in prev_nodes:
        line = str(index)+' '+str(prev_nodes[index]['type'])+' '+str(prev_nodes[index]['y'])+' '+\
               str(prev_nodes[index]['x'])+' '+str(prev_nodes[index]['z'])+' '+str(prev_nodes[index]['r'])+' '+\
               str(prev_nodes[index]['parent'])+'\n'
        out_f.write(line)
    out_f.close()


def write_to_pkl(fname, data):
    """
    HocCell objects maintain a nested dictionary specifying membrane mechanism parameters for each subcellular
    compartment. This method is used to save that dictionary to a .pkl file that can be read in during model
    specification or after parameter optimization.
    :param fname: str
    :param data: picklable object
    """
    output = open(fname, 'wb')
    pickle.dump(data, output, 2)
    output.close()


def read_from_pkl(fname):
    """
    HocCell objects maintain a nested dictionary specifying membrane mechanism parameters for each subcellular
    compartment. This method is used to load that dictionary from a .pkl file during model specification.
    :param fname: str
    :return: unpickled object
    """
    if os.path.isfile(fname):
        pkl_file = open(fname, 'rb')
        data = pickle.load(pkl_file)
        pkl_file.close()
        return data
    else:
        raise Exception('File: {} does not exist.'.format(fname))


def write_to_yaml(file_path, dict):
    """

    :param file_path: str (should end in '.yaml')
    :param dict: dict
    :return:
    """
    import yaml
    with open(file_path, 'w') as outfile:
        yaml.dump(dict, outfile, default_flow_style=False)


def read_from_yaml(file_path):
    """

    :param file_path: str (should end in '.yaml')
    :return:
    """
    import yaml
    if os.path.isfile(file_path):
        with open(file_path, 'r') as stream:
            data = yaml.load(stream)
        return data
    else:
        raise Exception('File: {} does not exist.'.format(file_path))


def combine_output_files(rec_file_list, new_rec_filename=None, local_data_dir=data_dir):
    """
    List contains names of files generated by "embarassingly parallel" execution of simulations on separate cores.
    This function combines the contents of the files into one .hdf5 file.
    :param rec_file_list: list
    :param new_rec_filename: str or None
    :param local_data_dir: str
    """
    if new_rec_filename is None:
        new_rec_filename = 'combined_output_'+datetime.datetime.today().strftime('%m%d%Y%H%M')
    new_f = h5py.File(local_data_dir+new_rec_filename+'.hdf5', 'w')
    simiter = 0
    for rec_filename in rec_file_list:
        old_f = h5py.File(local_data_dir+rec_filename+'.hdf5', 'r')
        for old_group in old_f.itervalues():
            new_f.copy(old_group, new_f, name=str(simiter))
            simiter += 1
        old_f.close()
    new_f.close()
    print 'Combined data in list of files and exported to: '+new_rec_filename+'.hdf5'
    return new_rec_filename


def combine_hdf5_file_paths(file_path_list, new_file_path=None):
    """
    List contains names of files generated by "embarassingly parallel" execution of simulations on separate cores.
    This function combines the contents of the files into one .hdf5 file.
    :param file_path_list: list of str (paths)
    :param new_file_path: str (path)
    """
    if new_file_path is None:
        raise ValueError('combine_output_file_paths: invalid file path provided: %s' % new_file_path)
    new_f = h5py.File(new_file_path, 'w')
    iter = 0
    for old_file_path in file_path_list:
        old_f = h5py.File(old_file_path, 'r')
        for old_group in old_f.itervalues():
            new_f.copy(old_group, new_f, name=str(iter))
            iter += 1
        old_f.close()
    new_f.close()
    print 'combine_output_file_paths: exported to file path: %s' % new_file_path


def time2index(tvec, start, stop):
    """
    When using adaptive time step (cvode), indices corresponding to specific time points cannot be calculated from a
    fixed dt. This method returns the indices closest to the duration bounded by the specified time points.
    :param tvec: :class:'numpy.array'
    :param start: float
    :param stop: float
    :return: tuple of int
    """
    left = np.where(tvec >= start)[0]
    if np.any(left):  # at least one value was found
        left = left[0]
    else:
        right = len(tvec) - 1  # just take the last two indices
        left = right - 1
        return left, right
    if tvec[left] >= stop:
        right = left
        left -= 1
        return left, right
    right = np.where(tvec <= stop)[0][-1]
    if right == left:
        left -= 1
    return left, right


def clean_axes(axes):
    """
    Remove top and right axes from pyplot axes object.
    :param axes:
    """
    if not type(axes) in [np.ndarray, list]:
        axes = [axes]
    elif type(axes) == np.ndarray:
        axes = axes.flatten()
    for axis in axes:
        axis.tick_params(direction='out')
        axis.spines['top'].set_visible(False)
        axis.spines['right'].set_visible(False)
        axis.get_xaxis().tick_bottom()
        axis.get_yaxis().tick_left()


def sort_str_list(str_list, seperator='_', end=None):
    """
    Given a list of filenames ending with (separator)int, sort the strings by increasing value of int.
    If there is a suffix at the end of the filename, provide it so it can be ignored.
    :param str_list: list of str
    :return: list of str
    """
    indexes = range(len(str_list))
    values = []
    for this_str in str_list:
        if end is not None:
            this_str = this_str.split(end)[0]
        this_value = int(this_str.split(seperator)[-1])
        values.append(this_value)
    indexes.sort(key=values.__getitem__)
    sorted_str_list = map(str_list.__getitem__, indexes)
    return sorted_str_list


def list_find(f, lst):
    index = None
    for i, x in enumerate(lst):
        if f(x):
            index = i
    return index


class Context(object):
    """
    A container replacement for global variables to be shared and modified by any function in a module.
    """
    def __init__(self, namespace_dict=None, **kwargs):
        self.update(namespace_dict, **kwargs)

    def update(self, namespace_dict=None, **kwargs):
        """
        Converts items in a dictionary (such as globals() or locals()) into context object internals.
        :param namespace_dict: dict
        """
        if namespace_dict is None:
            namespace_dict = {}
        namespace_dict.update(kwargs)
        for key, value in namespace_dict.iteritems():
            setattr(self, key, value)

    def __call__(self):
        return self.__dict__

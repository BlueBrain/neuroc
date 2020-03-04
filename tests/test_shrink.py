from pathlib import Path
import os

from morphio import SectionType
from morphio.mut import Morphology
from numpy.testing import (assert_array_almost_equal, assert_array_equal,
                           assert_equal)

from nose.tools import assert_dict_equal


from neuroc.axon_shrinker.shrink import *

DATA_PATH = Path(Path(__file__).parent, 'data')

def path(file):
    return str(Path(DATA_PATH, file))

def _neuron():
    return Morphology(path('simple.asc'))


def _get_axon(neuron):
    return next(sec for sec in neuron.root_sections if sec.type == SectionType.axon)


def test_section_path_length():
    neuron = _neuron()

    assert_dict_equal(section_path_lengths(neuron.root_sections[0]),
                       {0: 5.0, 1: 10.0, 2: 11.0})


    assert_dict_equal(section_path_lengths(neuron.root_sections[1]),
                      {3: 4.123105525970459, 4: 10.123105525970459, 5: 9.123105525970459})


def test_cut_and_graft_coordinate():
    neuron = _neuron()

    pos1 = {'y_min': 0, 'y_max': 10}
    pos2 = {'y_min': 100, 'y_max': 200}
    assert_array_equal(cut_and_graft_coordinate({'dendrite': pos1, 'axon': pos2, }),
                       [True, 10, 100])

    assert_array_equal(cut_and_graft_coordinate({'dendrite': pos2, 'axon': pos1}),
                       [False, 100, 10])


def test_cut_branch():
    neuron = _neuron()

    start_cut = _get_axon(neuron)
    cut_branch(neuron, upward=False, start_cut=start_cut, y_start_cut=-2.5)
    new_axon = next(sec for sec in neuron.root_sections if sec.type == SectionType.axon)
    assert_array_equal(new_axon.points,
                       [[0, 0, 1], [0, -2.5, 0.375]])
    assert_array_equal(new_axon.diameters,
                       [1, 1.625])


def test_add_vertical_segment():
    neuron = _neuron()

    end_axon = _get_axon(neuron).children[1]
    add_vertical_segment(neuron, end_axon, 11.5)
    assert_equal(len(end_axon.children), 1)

    assert_array_equal(end_axon.children[0].points,
                       [end_axon.points[-1],
                        end_axon.points[-1] + [0, 11.5, 0]])


def test_graft_branch():
    neuron = _neuron()
    dendrite = next(sec for sec in neuron.root_sections if sec.type != SectionType.axon)

    end_axon = _get_axon(neuron).children[1]

    y_diff = graft_branch(True, end_axon, dendrite, -4)

    assert_equal(y_diff, -4)
    assert_equal(len(end_axon.children), 1)
    grafted_dendrite = end_axon.children[0]
    assert_array_equal(grafted_dendrite.points,
                       [[-5, -4, 0], [-5, 1, 0]])

    assert_equal(len(grafted_dendrite.children), 2)


def test_cut_and_graft():
    neuron = _neuron()
    new_neuron, y_diff = cut_and_graft(path('neuron.asc'),
                                       upward=True,
                                       y_start_cut=1.5,
                                       y_start_graft=4.6,
                                       height=103)

    assert_equal(len(new_neuron.root_sections), 1)
    root = new_neuron.root_sections[0]
    assert_array_equal(root.points,
                       [[0, 0, 0], [0, 1, 0], [0, 1.5, 0]])

    assert_equal(len(root.children), 1)
    vertical = root.children[0]
    assert_equal(vertical.points,
                 [[0, 1.5, 0], [0, 104.5, 0]])

    assert_array_equal(len(vertical.children), 1)
    graft_root = vertical.children[0]
    assert_array_equal(graft_root.points,
                       np.array([[0, 104.5, 0],
                                 [0, 104.9, 0],
                                 [0, 105.4, 0]], dtype=np.float32))

    assert_equal(len(graft_root.children), 2)
    children = graft_root.children
    assert_array_equal(children[0].points,
                       np.array([[0, 105.4, 0], [0, 105.9, 1]],
                                dtype=np.float32))
    assert_array_equal(children[1].points,
                       np.array([[0, 105.4, 0], [0, 106.9, 2]],
                                dtype=np.float32))


def test_cut_and_graft_downward():
    neuron = _neuron()
    new_neuron, y_diff = cut_and_graft(path('neuron_downward.asc'),
                                       upward=False,
                                       y_start_cut=-1.5,
                                       y_start_graft=-4.6,
                                       height=-103)

    assert_equal(len(new_neuron.root_sections), 1)
    root = new_neuron.root_sections[0]
    assert_array_equal(root.points,
                       [[0, 0, 0], [0, -1, 0], [0, -1.5, 0]])

    assert_equal(len(root.children), 1)
    vertical = root.children[0]
    assert_equal(vertical.points,
                 [[0, -1.5, 0], [0, -104.5, 0]])

    assert_array_equal(len(vertical.children), 1)
    graft_root = vertical.children[0]
    assert_array_equal(graft_root.points,
                       np.array([[0, -104.5, 0],
                                 [0, -104.9, 0],
                                 [0, -105.4, 0]], dtype=np.float32))

    assert_equal(len(graft_root.children), 2)
    children = graft_root.children
    assert_array_equal(children[0].points,
                       np.array([[0, -105.4, 0], [0, -105.9, 1]],
                                dtype=np.float32))
    assert_array_equal(children[1].points,
                       np.array([[0, -105.4, 0], [0, -106.9, 2]],
                                dtype=np.float32))

def test_cut_axon_end():
    neuron = cut_axon_end(path('axon.asc'), -9)
    axon = neuron.root_sections[0]
    assert_array_equal(axon.points,
                       np.array([[ 0.,  0.,  1.],
                                 [ 0., -4.,  0.]], dtype=np.float32))

    assert_array_equal(axon.children[0].points,
                       np.array([[ 0., -4.,  0.],
                                 [6, -6, 0]], dtype=np.float32))

    assert_array_equal(axon.children[1].points,
                       np.array([[ 0., -4.,  0.],
                                 [-5, -6, 0],
                                 [-5, -8, 0],
                                 [-5, -9, 0]], dtype=np.float32))

def test_cut_axon_end_backward():
    cut_axon_end(path('axon-to-cut-end.h5'), 500)

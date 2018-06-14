import os

import morphio
from morphio.mut import Morphology
from numpy.testing import (assert_array_almost_equal, assert_array_equal,
                           assert_equal)

from neuroc.axon_shrinker.shrink import *

_path = os.path.dirname(os.path.abspath(__file__))


def _neuron():
    return Morphology(os.path.join(_path, 'simple.asc'))


def _get_axon(neuron):
    return next(sec for sec in neuron.root_sections if sec.type == morphio.axon)


def test_section_path_length():
    neuron = _neuron()

    assert_array_equal(section_path_lengths(neuron, neuron.root_sections[0]).values(),
                       [5, 10, 11])
    assert_array_almost_equal(section_path_lengths(neuron, neuron.root_sections[1]).values(),
                              [4.123106, 10.123106, 9.123106])


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
    new_axon = next(sec for sec in neuron.root_sections if sec.type == morphio.axon)
    assert_array_equal(new_axon.points,
                       [[0, 0, 1], [0, -2.5, 0.375]])
    assert_array_equal(new_axon.diameters,
                       [1, 1.625])


def test_add_vertical_segment():
    neuron = _neuron()

    end_axon = neuron.children(_get_axon(neuron))[1]
    add_vertical_segment(neuron, end_axon, 11.5)
    assert_equal(len(neuron.children(end_axon)), 1)

    assert_array_equal(neuron.children(end_axon)[0].points,
                       [end_axon.points[-1],
                        end_axon.points[-1] + [0, 11.5, 0]])


def test_graft_branch():
    neuron = _neuron()
    dendrite = next(sec for sec in neuron.root_sections if sec.type != morphio.axon)

    end_axon = neuron.children(_get_axon(neuron))[1]

    y_diff = graft_branch(neuron, neuron, True, end_axon, dendrite, -4)

    assert_equal(y_diff, -4)
    assert_equal(len(neuron.children(end_axon)), 1)
    grafted_dendrite = neuron.children(end_axon)[0]
    assert_array_equal(grafted_dendrite.points,
                       [[-5, -4, 0], [-5, 1, 0]])

    assert_equal(len(neuron.children(grafted_dendrite)), 2)


def test_cut_and_graft():
    neuron = _neuron()
    new_neuron, y_diff = cut_and_graft(os.path.join(_path, 'neuron.asc'),
                                       upward=True,
                                       y_start_cut=1.5,
                                       y_start_graft=4.6,
                                       height=103)

    assert_equal(len(new_neuron.root_sections), 1)
    root = new_neuron.root_sections[0]
    assert_array_equal(root.points,
                       [[0, 0, 0], [0, 1, 0], [0, 1.5, 0]])

    assert_equal(len(new_neuron.children(root)), 1)
    vertical = new_neuron.children(root)[0]
    assert_equal(vertical.points,
                 [[0, 1.5, 0], [0, 104.5, 0]])

    assert_array_equal(len(new_neuron.children(vertical)), 1)
    graft_root = new_neuron.children(vertical)[0]
    assert_array_equal(graft_root.points,
                       np.array([[0, 104.5, 0],
                                 [0, 104.9, 0],
                                 [0, 105.4, 0]], dtype=np.float32))

    assert_equal(len(new_neuron.children(graft_root)), 2)
    children = new_neuron.children(graft_root)
    assert_array_equal(children[0].points,
                       np.array([[0, 105.4, 0], [0, 105.9, 1]],
                                dtype=np.float32))
    assert_array_equal(children[1].points,
                       np.array([[0, 105.4, 0], [0, 106.9, 2]],
                                dtype=np.float32))


def test_cut_and_graft_downward():
    neuron = _neuron()
    new_neuron, y_diff = cut_and_graft(os.path.join(_path, 'neuron_downward.asc'),
                                       upward=False,
                                       y_start_cut=-1.5,
                                       y_start_graft=-4.6,
                                       height=-103)

    assert_equal(len(new_neuron.root_sections), 1)
    root = new_neuron.root_sections[0]
    assert_array_equal(root.points,
                       [[0, 0, 0], [0, -1, 0], [0, -1.5, 0]])

    assert_equal(len(new_neuron.children(root)), 1)
    vertical = new_neuron.children(root)[0]
    assert_equal(vertical.points,
                 [[0, -1.5, 0], [0, -104.5, 0]])

    assert_array_equal(len(new_neuron.children(vertical)), 1)
    graft_root = new_neuron.children(vertical)[0]
    assert_array_equal(graft_root.points,
                       np.array([[0, -104.5, 0],
                                 [0, -104.9, 0],
                                 [0, -105.4, 0]], dtype=np.float32))

    assert_equal(len(new_neuron.children(graft_root)), 2)
    children = new_neuron.children(graft_root)
    assert_array_equal(children[0].points,
                       np.array([[0, -105.4, 0], [0, -105.9, 1]],
                                dtype=np.float32))
    assert_array_equal(children[1].points,
                       np.array([[0, -105.4, 0], [0, -106.9, 2]],
                                dtype=np.float32))

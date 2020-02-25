"""Sample nosetest file."""
import warnings
from pathlib import Path
from tempfile import TemporaryDirectory

from neurom import load_neuron, COLS
from nose.tools import assert_dict_equal, assert_equal, assert_raises
from numpy.testing import assert_almost_equal, assert_array_almost_equal

import neuroc.rat_to_human as test_module
from morphio.mut import Morphology

DATA = Path(__file__).parent / 'data'
RAT_PATH = DATA / 'rat-cells'
HUMAN_PATH = DATA / 'human-cells'
INH_MAPPING_PATH = DATA / 'inh-mapping.yaml'
EXC_MAPPING_PATH = DATA / 'exc-mapping.yaml'

EXCITATORY_PATH = HUMAN_PATH / 'Excitatory/'
INHIBITORY_PATH = HUMAN_PATH / 'Inhibitory/'


def test__find_filepath():
    assert_equal(test_module._find_filepath(DATA, 'Neuron'), Path(DATA, 'Neuron.swc'))
    assert_equal(test_module._find_filepath(DATA, 'miles.davis'), None)


def test_human_cells_per_layer_mtype():
    assert_equal(set(test_module.human_cells_per_layer_mtype(EXCITATORY_PATH).keys()),
                 {'L1', 'L2', 'L3', 'L4', 'L5', 'L6'})


def test_validate_folders():
    wrong_path = DATA
    assert_raises(ValueError, test_module.validate_folders, wrong_path, RAT_PATH)
    assert_raises(ValueError, test_module.validate_folders, EXCITATORY_PATH, wrong_path)


def test_morphology_filenames():
    expected = {'L1_DAC': [RAT_PATH / 'neuron1.swc'],
                'L4_TPC': [RAT_PATH / 'neuron2.swc'],
                'some-mtype': [RAT_PATH / 'neuron3.swc']}
    assert_dict_equal(test_module.morphology_filenames(RAT_PATH), expected)

def test_extensions():
    assert_almost_equal(test_module.dendritic_y_std(load_neuron(RAT_PATH / 'neuron1.swc')),
                        2.165063509461097)
    assert_almost_equal(test_module.dendritice_radial_std(load_neuron(RAT_PATH / 'neuron1.swc')),
                        2.7726341266023544)
    assert_almost_equal(test_module.dendritic_diameter(load_neuron(RAT_PATH / 'neuron1.swc')), 2.5)

    assert_almost_equal(test_module.dendritic_y_std(load_neuron(DATA / 'Neuron.swc')),
                        28.280295211604614)
    assert_almost_equal(test_module.dendritice_radial_std(load_neuron(DATA / 'Neuron.swc')),
                        22.01310793432494)
    assert_almost_equal(test_module.dendritic_diameter(load_neuron(DATA / 'Neuron.swc')),
                        1.2016682319539242)

def test_scaling_factors():

    group1 = [RAT_PATH / 'neuron1.swc', RAT_PATH / 'neuron1.swc', RAT_PATH / 'neuron1.swc']
    group2 = [DATA / 'Neuron.swc', DATA / 'Neuron.swc', DATA / 'Neuron.swc']
    assert_array_almost_equal(test_module.scaling_factors(group1, group2,
                                                          [test_module.dendritice_radial_std,
                                                           test_module.dendritic_y_std,
                                                           test_module.dendritic_diameter]),
                              [0.125954, 0.076557, 2.080441])


def test_scale_diameter():
    neuron = Morphology(RAT_PATH / 'neuron1.swc')
    test_module.scale_diameter(neuron, 2.)
    assert_array_almost_equal(neuron.section(0).diameters,
                              [4, 4])


def test_morphologies_in_layer():
    assert test_module.morphologies_in_layer(RAT_PATH, 'L3') == [
        RAT_PATH / 'neuron3.swc'
    ]

def test_mtype_matcher():
    with warnings.catch_warnings(record=True):
        groups = test_module.mtype_matcher(INHIBITORY_PATH, RAT_PATH, INH_MAPPING_PATH)

    def assert_same_groups(actual, expected):
        assert len(actual) == len(expected)
        for (act_left, act_right), (exp_left, exp_right) in zip(actual, expected):
            assert {str(path) for path in act_left} == {str(path) for path in exp_left}
            assert_equal({str(path) for path in act_right},
                         {str(path) for path in exp_right})

    assert_same_groups(groups,
                       [([INHIBITORY_PATH / 'L1/PSC_some-cell-name.swc'], [RAT_PATH / 'neuron1.swc']),
                        ([INHIBITORY_PATH / 'L2/AC_some-cell-name.swc'], []),

                      # The following line is the result of the "ALL" matching
                      ([INHIBITORY_PATH / 'L3/AC_some-cell-name.swc'], [RAT_PATH / 'neuron3.swc']),
                        ([INHIBITORY_PATH / 'L4/BTC_some-cell-name.swc'], []),
                        ([INHIBITORY_PATH / 'L5/MPC_some-cell-name.swc'], [])])


    with warnings.catch_warnings(record=True):
        groups = test_module.mtype_matcher(EXCITATORY_PATH, RAT_PATH, EXC_MAPPING_PATH)

    assert_same_groups(groups,
                       [([EXCITATORY_PATH / 'L1/PSC_some-cell-name.swc'], []),
                        ([EXCITATORY_PATH / 'L2/AC_some-cell-name.swc'], []),
                        ([EXCITATORY_PATH / 'L3/AC_some-cell-name.swc'], []),
                        ([EXCITATORY_PATH / 'L4/BTC_some-cell-name.swc'], [RAT_PATH / 'neuron2.swc']),
                        ([EXCITATORY_PATH / 'L5/MPC_some-cell-name.swc'], [])])

    assert_raises(ValueError, test_module.mtype_matcher, INHIBITORY_PATH, RAT_PATH, DATA / 'broken-mapping.yaml')

def test_scale_all_cells():
    with TemporaryDirectory('test-scale-rat-cells') as output_folder:
        output_folder = Path(output_folder)
        with warnings.catch_warnings(record=True):
            test_module.scale_all_cells(INHIBITORY_PATH, RAT_PATH, INH_MAPPING_PATH, output_folder)
        assert_equal(list(output_folder.rglob('*')),
                     [
                         output_folder / 'neuron1_-_Y-Scale_2.0_-_XZ-Scale_2.0_-_Diam-Scale_3.0.h5',
                         output_folder / 'neuron3_-_Y-Scale_1.0_-_XZ-Scale_1.0_-_Diam-Scale_1.0.h5',
                     ])


def test_scale_single_coordinates():
    orig_neuron = Morphology(DATA / 'Neuron.swc')
    neuron = Morphology(orig_neuron)
    scaling_value = 4
    test_module.scale_coordinates(neuron, scaling_value, COLS.Y)
    points, orig_points = neuron.section(0).points, orig_neuron.section(0).points
    assert_array_almost_equal(points[:, COLS.Y], orig_points[:, COLS.Y] * scaling_value)
    assert_array_almost_equal(points[:, COLS.XZ], orig_points[:, COLS.XZ])


def test_scale_double_coordinates():
    orig_neuron = Morphology(DATA / 'Neuron.swc')
    neuron = Morphology(orig_neuron)
    scaling_value = 4
    test_module.scale_coordinates(neuron, scaling_value, COLS.XZ)
    points, orig_points = neuron.section(0).points, orig_neuron.section(0).points
    assert_array_almost_equal(points[:, COLS.XZ], orig_points[:, COLS.XZ] * scaling_value)
    assert_array_almost_equal(points[:, COLS.Y], orig_points[:, COLS.Y])


def test_scale_one_cells():
    filename = DATA / 'rp110711_C3_idA.h5'
    orig = Morphology(filename)
    neuron = test_module.scale_one_cell(filename, 2.5, 1, 1)

    s1, s2 = orig.section(0), neuron.section(0)
    assert_array_almost_equal(s2.points[:, COLS.Y] / s1.points[:, COLS.Y], 2.5)
    assert_array_almost_equal(s2.points[:, COLS.XZ] / s1.points[:, COLS.XZ], 1)
    assert_array_almost_equal(s2.diameters / s1.diameters, 1)

    orig = Morphology(filename)
    neuron = test_module.scale_one_cell(filename, 3.3, 4.5, 12.7)

    s1, s2 = orig.section(0), neuron.section(0)
    assert_array_almost_equal(s2.points[:, COLS.Y] / s1.points[:, COLS.Y], 3.3)
    assert_array_almost_equal(s2.points[:, COLS.XZ] / s1.points[:, COLS.XZ], 4.5)
    assert_array_almost_equal(s2.diameters / s1.diameters, 12.7)

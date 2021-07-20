"""Sample nosetest file."""
import warnings
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd
from morphio.mut import Morphology
from neurom import COLS, load_neuron
import pytest
from numpy.testing import (assert_almost_equal, assert_array_almost_equal,
                           assert_array_equal)
from pandas.testing import assert_frame_equal

import neuroc.rat_to_human as tested

DATA = Path(__file__).parent / 'data'
RAT_PATH = DATA / 'rat-cells/h5'
HUMAN_PATH = DATA / 'human-cells/h5'
MAPPING_PATH = DATA / 'mapping.yaml'


def test__find_filepath():
    assert tested._find_filepath(DATA, 'Neuron') == Path(DATA, 'Neuron.swc')
    with pytest.raises(ValueError):
        tested._find_filepath(DATA, 'miles.davis')


def test_extensions():
    assert_almost_equal(tested.dendritic_y_std(load_neuron(RAT_PATH / 'neuron1.swc')),
                        2.165063509461097)
    assert_almost_equal(tested.dendritice_radial_std(load_neuron(RAT_PATH / 'neuron1.swc')),
                        2.7726341266023544)
    assert_almost_equal(tested.dendritic_diameter(load_neuron(RAT_PATH / 'neuron1.swc')), 2.5)

    assert_almost_equal(tested.dendritic_y_std(load_neuron(DATA / 'Neuron.swc')),
                        28.28029441833496)
    assert_almost_equal(tested.dendritice_radial_std(load_neuron(DATA / 'Neuron.swc')),
                        22.013107299804688)
    assert_almost_equal(tested.dendritic_diameter(load_neuron(DATA / 'Neuron.swc')),
                        1.2016682624816895)


def test_scaling_factors():

    group1 = [RAT_PATH / 'neuron1.swc', RAT_PATH / 'neuron1.swc', RAT_PATH / 'neuron1.swc']
    group2 = [DATA / 'Neuron.swc', DATA / 'Neuron.swc', DATA / 'Neuron.swc']
    assert_array_almost_equal(tested.scaling_factors(group1, group2,
                                                     [tested.dendritice_radial_std,
                                                      tested.dendritic_y_std,
                                                      tested.dendritic_diameter]),
                              [0.125954, 0.076557, 2.080441])


def test_scale_diameter():
    neuron = Morphology(RAT_PATH / 'neuron1.swc')
    tested.scale_diameter(neuron, 2.)
    assert_array_almost_equal(neuron.section(0).diameters,
                              [4, 4])


def test_mtype_matcher():
    with warnings.catch_warnings(record=True):
        groups = tested.mtype_matcher(HUMAN_PATH / '../neurondb.dat',
                                      RAT_PATH / '../neurondb.dat',
                                      MAPPING_PATH)

    human, rat = next(groups)
    assert_array_equal(human.name, ['AC_some-cell-name'])
    assert_array_equal(rat.name, ['neuron1'])

    human, rat = next(groups)
    assert_array_equal(human.name, ['PSC_some-cell-name'])
    assert_array_equal(rat.name, ['neuron2', 'neuron3', 'neuron4'])


def test_scale_single_coordinates():
    orig_neuron = Morphology(DATA / 'Neuron.swc')
    neuron = Morphology(orig_neuron)
    scaling_value = 4
    tested.scale_coordinates(neuron, scaling_value, COLS.Y)
    points, orig_points = neuron.section(0).points, orig_neuron.section(0).points
    assert_array_almost_equal(points[:, COLS.Y], orig_points[:, COLS.Y] * scaling_value)
    assert_array_almost_equal(points[:, COLS.XZ], orig_points[:, COLS.XZ])


def test_scale_double_coordinates():
    orig_neuron = Morphology(DATA / 'Neuron.swc')
    neuron = Morphology(orig_neuron)
    scaling_value = 4
    tested.scale_coordinates(neuron, scaling_value, COLS.XZ)
    points, orig_points = neuron.section(0).points, orig_neuron.section(0).points
    assert_array_almost_equal(points[:, COLS.XZ], orig_points[:, COLS.XZ] * scaling_value)
    assert_array_almost_equal(points[:, COLS.Y], orig_points[:, COLS.Y])


def test_scale_one_cells():
    filename = DATA / 'rp110711_C3_idA.h5'
    orig = Morphology(filename)
    neuron = tested.scale_one_cell(filename, 2.5, 1, 1)

    s1, s2 = orig.section(0), neuron.section(0)
    assert_array_almost_equal(s2.points[:, COLS.Y] / s1.points[:, COLS.Y], 2.5)
    assert_array_almost_equal(s2.points[:, COLS.XZ] / s1.points[:, COLS.XZ], 1)
    assert_array_almost_equal(s2.diameters / s1.diameters, 1)

    orig = Morphology(filename)
    neuron = tested.scale_one_cell(filename, 3.3, 4.5, 12.7)

    s1, s2 = orig.section(0), neuron.section(0)
    assert_array_almost_equal(s2.points[:, COLS.Y] / s1.points[:, COLS.Y], 3.3)
    assert_array_almost_equal(s2.points[:, COLS.XZ] / s1.points[:, COLS.XZ], 4.5)
    assert_array_almost_equal(s2.diameters / s1.diameters, 12.7)


def test_scale_all_cells():
    with TemporaryDirectory('test-scale-rat-cells') as output_folder:
        output_folder = Path(output_folder)
        with warnings.catch_warnings(record=True):
            tested.scale_all_cells(HUMAN_PATH / '../neurondb.dat',
                                   RAT_PATH / '../neurondb.dat',
                                   MAPPING_PATH, output_folder)
        assert (set(output_folder.rglob('*')) ==
                     {
                         Path(
                             output_folder, 'neuron3_-_Y-Scale_2.0_-_XZ-Scale_2.0_-_Diam-Scale_3.0.h5'),
                         Path(
                             output_folder, 'neuron2_-_Y-Scale_2.0_-_XZ-Scale_2.0_-_Diam-Scale_3.0.h5'),
                         Path(output_folder, 'neurondb.csv'),
                         Path(
                             output_folder, 'neuron4_-_Y-Scale_2.0_-_XZ-Scale_2.0_-_Diam-Scale_3.0.h5'),
                         Path(output_folder, 'neuron1_-_Y-Scale_1.0_-_XZ-Scale_1.0_-_Diam-Scale_1.0.h5')
        })

        df = pd.read_csv(output_folder / 'neurondb.csv', index_col=False)
        expected = pd.read_csv(DATA / 'expected-metadata.csv')
        assert_frame_equal(df, expected)


def test_extended_neurondb():
    assert_frame_equal(tested.extended_neurondb(RAT_PATH / '../neurondb.dat'),
                       pd.DataFrame({'name': ['neuron1', 'neuron2', 'neuron3', 'neuron4'],
                                     'layer': ['L1', 'L5', 'L5', 'L5'],
                                     'mtype': ['L1_DAC', 'L5_PC', 'L5_PC', 'L5_HAC'],
                                     'path': [Path(RAT_PATH, 'neuron1.swc'),
                                              Path(RAT_PATH, 'neuron2.swc'),
                                              Path(RAT_PATH, 'neuron3.swc'),
                                              Path(RAT_PATH, 'neuron4.swc')]}))

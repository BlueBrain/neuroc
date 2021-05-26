from pathlib import Path

import numpy as np
from morph_tool import diff
from morphio import Morphology as ImmutMorphology
from morphio.mut import Morphology
from nose.tools import ok_
from numpy.testing import (assert_array_almost_equal, assert_array_equal,
                           assert_equal)

import neuroc.scale as tested
from neuroc.scale import RotationParameters, ScaleParameters

DATA_PATH = Path(Path(__file__).parent, 'data')


def path(file):
    return str(Path(DATA_PATH, file))


SIMPLE_PATH = path('simple.asc')
NEURON_PATH = path('Neuron.swc')


def test_segment_vector():
    neuron = Morphology(SIMPLE_PATH)
    assert_array_equal(tested._segment_vectors(neuron.section(0), prepend_null_vector=False),
                       [[0, 5, 0]])

    assert_array_equal(tested._segment_vectors(neuron.section(0), prepend_null_vector=True),
                       [[0, 0, 0], [0, 5, 0]])


def test_rotational_jitter():
    neuron = Morphology(NEURON_PATH)
    tested.rotational_jitter(neuron, RotationParameters(mean_angle=0, std_angle=0, numberpoint=5))
    ok_(not diff(NEURON_PATH, neuron))

    neuron = Morphology(NEURON_PATH)
    tested.rotational_jitter(neuron, RotationParameters(mean_angle=360, std_angle=0, numberpoint=5))
    ok_(not diff(NEURON_PATH, neuron))

    neuron = Morphology(SIMPLE_PATH)
    tested.rotational_jitter(neuron, RotationParameters(numberpoint=5, mean_angle=90., std_angle=0))

    # The parent section is oriented along Y so this is a rotation around Y
    # For reference, the original section is: [[ 0., 5., 0.], [-5., 5., 0.]]
    assert_array_almost_equal(neuron.section(1).points, [[0, 5, 0], [0, 5, 5]])

    np.random.seed(0)
    neuron = Morphology(SIMPLE_PATH)
    tested.rotational_jitter(neuron, RotationParameters(numberpoint=5, mean_angle=90., std_angle=0.1))
    assert_array_almost_equal(neuron.section(1).points, [[0, 5, 0], [0.015394, 5, 4.999976]])

    rng = np.random.default_rng(0)
    neuron = Morphology(SIMPLE_PATH)
    tested.rotational_jitter(
        neuron,
        RotationParameters(numberpoint=5, mean_angle=90., std_angle=0.1),
        rng=rng,
    )
    # Note: the random numbers generated after np.random.seed(0) and np.random.default_rng(0)
    # are different so it is expected that the arrays are slightly different.
    assert_array_almost_equal(neuron.section(1).points, [[0, 5, 0], [0.016297, 5, 4.999973]])


def test_no_scaling():
    neuron = Morphology(SIMPLE_PATH)
    tested.scale_morphology(neuron, ScaleParameters(), ScaleParameters())
    ok_(not diff(SIMPLE_PATH, neuron))


def test_section_scaling():
    # Scale by 100%
    neuron = Morphology(SIMPLE_PATH)
    section = neuron.section(0)
    tested.scale_section(section, ScaleParameters(mean=2))

    # section(0) scaled
    assert_array_almost_equal(neuron.section(0).points,
                              [[0, 0, 0], [0, 10, 0]])
    # section(1) tranlated but not scaled
    assert_array_almost_equal(neuron.section(1).points,
                              [[0, 10, 0], [-5, 10, 0]])


def test_morphology_scaling_section_param_only():
    # Scale by 100%
    neuron = Morphology(SIMPLE_PATH)
    tested.scale_morphology(neuron, ScaleParameters(mean=2))
    assert_array_almost_equal(neuron.section(0).points,
                              [[0, 0, 0], [0, 10, 0]])
    assert_array_almost_equal(neuron.section(1).points,
                              [[0, 10, 0], [-10, 10, 0]])

    # Scaling only on X axis
    neuron = Morphology(SIMPLE_PATH)
    tested.scale_morphology(neuron, ScaleParameters(mean=2, axis=0))
    assert_array_almost_equal(neuron.section(0).points,
                              [[0, 0, 0], [0, 5, 0]])
    assert_array_almost_equal(neuron.section(1).points,
                              [[0, 5, 0], [-10, 5., 0]])

    # Attempt at scaling by -200% but minimum scaling factor is 1% of original length
    neuron = Morphology(SIMPLE_PATH)
    tested.scale_morphology(neuron, ScaleParameters(mean=-2))
    assert_array_almost_equal(neuron.section(0).points,
                              [[0, 0, 0], [0, 0.05, 0]])
    assert_array_almost_equal(neuron.section(1).points,
                              [[0, 0.05, 0], [-0.05, 0.05, 0]])


def test_morphology_scaling_segment_param_only():
    neuron = Morphology(SIMPLE_PATH)
    tested.scale_morphology(neuron, ScaleParameters(), ScaleParameters(mean=2))
    assert_array_almost_equal(neuron.section(0).points,
                              [[0, 0, 0], [0, 10, 0]])
    assert_array_almost_equal(neuron.section(1).points,
                              [[0, 10, 0], [-10, 10, 0]])

    neuron = Morphology(SIMPLE_PATH)
    np.random.seed(0)
    tested.scale_morphology(neuron, ScaleParameters(), ScaleParameters(mean=2, std=0.5))
    assert_array_almost_equal(neuron.section(0).points,
                              [[0, 0, 0], [0., 9.621607, 0.]])
    assert_array_almost_equal(neuron.section(1).points,
                              [[0., 9.621607, 0.], [-11.000393, 9.621607, 0.]])

    neuron = Morphology(SIMPLE_PATH)
    rng = np.random.default_rng(0)
    tested.scale_morphology(neuron, ScaleParameters(), ScaleParameters(mean=2, std=0.5), rng=rng)
    assert_array_almost_equal(neuron.section(0).points,
                              [[0, 0, 0], [0., 9.669738, 0.]])
    assert_array_almost_equal(neuron.section(1).points,
                              [[0., 9.669738, 0.], [-7.864761, 9.669738, 0.]])

    # Scaling only on Y axis
    neuron = Morphology(SIMPLE_PATH)
    tested.scale_morphology(neuron, ScaleParameters(), ScaleParameters(mean=2, axis=1))
    assert_array_almost_equal(neuron.section(0).points,
                              [[0, 0, 0], [0, 10, 0]])
    assert_array_almost_equal(neuron.section(1).points,
                              [[0, 10, 0], [-5, 10, 0]])


def test_principal_direction():
    neuron = Morphology(SIMPLE_PATH)
    assert_array_almost_equal(tested._principal_direction(
        neuron.section(0)), [0.998492, -0.0549, 0.])

    neuron.section(1).points = [[0, 0, 0], [1, 1, 0], [-1, 1, 0]]
    assert_array_equal(tested._principal_direction(neuron.section(1)), [1, 0, 0])

    neuron = Morphology(path('simple-with-duplicate.asc'))
    assert_array_almost_equal(tested._principal_direction(neuron.section(1)),
                              [0.098538, 0.995133, 0.      ])


def test_iter_clones():
    clone = next(tested.yield_clones(NEURON_PATH, RotationParameters(30, 0, 5)))
    expected = ImmutMorphology(path('neuron_rotation_30_degree.swc'))
    ok_(not diff(clone, expected, atol=1e-3))

    # Test with segment_scaling
    it = tested.yield_clones(SIMPLE_PATH, segment_scaling=ScaleParameters(mean=2, std=1), seed=0)
    clone = next(it)
    assert_array_almost_equal(clone.section(0).points,
                              [[0, 0, 0], [0., 5.113611, 0.]])
    assert_array_almost_equal(clone.section(1).points,
                              [[0., 5.113611, 0.], [-12.000786, 5.113611, 0.]])

    clone = next(it)
    assert_array_almost_equal(clone.section(0).points,
                              [[0., 0., 0.],
                               [0., 8.974209, 0.]])
    assert_array_almost_equal(clone.section(1).points,
                              [[0., 8.974209, 0.],
                               [-10.608376, 8.974209, 0.]])

    # Test with section_scaling
    it = tested.yield_clones(SIMPLE_PATH, section_scaling=ScaleParameters(mean=2))
    clone = next(it)
    assert_array_almost_equal(clone.section(0).points,
                              [[0, 0, 0], [0, 10, 0]])
    assert_array_almost_equal(clone.section(1).points,
                              [[0, 10, 0], [-10, 10, 0]])

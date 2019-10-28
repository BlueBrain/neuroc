from pathlib import Path
from tempfile import TemporaryDirectory
import os

from morph_tool import diff
from morphio.mut import Morphology
from morphio import Morphology as ImmutMorphology
import numpy as np
from numpy.testing import (assert_array_almost_equal, assert_array_equal,
                           assert_equal)
from nose.tools import ok_

from neuroc.jitter import (rotational_jitter, RotationParameters, ScaleParameters,
                           scale_morphology, create_clones, _principal_direction, _segment_vectors)

DATA_PATH = Path(Path(__file__).parent, 'data')

def path(file):
    return str(Path(DATA_PATH, file))

SIMPLE_PATH = path('simple.asc')
NEURON_PATH = path('Neuron.swc')

def test_segment_vector():
    neuron = Morphology(SIMPLE_PATH)
    assert_array_equal(_segment_vectors(neuron.section(0), prepend_null_vector=False),
                       [[0, 5, 0]])

    assert_array_equal(_segment_vectors(neuron.section(0), prepend_null_vector=True),
                       [[0, 0, 0], [0, 5, 0]])

def test_rotational_jitter():
    neuron = Morphology(NEURON_PATH)
    rotational_jitter(neuron, RotationParameters(mean_angle=0, std_angle=0, numberpoint=5))
    ok_(not diff(NEURON_PATH, neuron))

    neuron = Morphology(NEURON_PATH)
    rotational_jitter(neuron, RotationParameters(mean_angle=360, std_angle=0, numberpoint=5))
    ok_(not diff(NEURON_PATH, neuron))

    neuron = Morphology(SIMPLE_PATH)
    rotational_jitter(neuron, RotationParameters(numberpoint=5, mean_angle=90., std_angle=0))

    # The parent section is oriented along Y so this is a rotation around Y
    # For reference, the original section is: [[ 0., 5., 0.], [-5., 5., 0.]]
    assert_array_almost_equal(neuron.section(1).points, [[0, 5, 0], [0, 5, 5]])


def test_no_scaling():
    neuron = Morphology(SIMPLE_PATH)
    scale_morphology(neuron, ScaleParameters(), ScaleParameters())
    ok_(not diff(SIMPLE_PATH, neuron))

def test_section_scaling():
    # Scale by 100%
    neuron = Morphology(SIMPLE_PATH)
    scale_morphology(neuron, ScaleParameters(), ScaleParameters(mean=2))
    assert_array_almost_equal(neuron.section(0).points,
                              [[0, 0, 0], [0, 10, 0]])
    assert_array_almost_equal(neuron.section(1).points,
                              [[0, 10, 0], [-10, 10, 0]])

    # Scaling only on X axis
    neuron = Morphology(SIMPLE_PATH)
    scale_morphology(neuron, ScaleParameters(), ScaleParameters(mean=2, axis=0))
    assert_array_almost_equal(neuron.section(0).points,
                              [[0, 0, 0], [0, 5, 0]])
    assert_array_almost_equal(neuron.section(1).points,
                              [[0, 5, 0], [-10, 5., 0]])

    # Attempt at scaling by -200% but minimum scaling factor is 1% of original length
    neuron = Morphology(SIMPLE_PATH)
    scale_morphology(neuron, ScaleParameters(), ScaleParameters(mean=-2))
    assert_array_almost_equal(neuron.section(0).points,
                              [[0, 0, 0], [0, 0.05, 0]])
    assert_array_almost_equal(neuron.section(1).points,
                              [[0, 0.05, 0], [-0.05, 0.05, 0]])

def test_segment_scaling():
    neuron = Morphology(SIMPLE_PATH)
    scale_morphology(neuron, ScaleParameters(mean=2), ScaleParameters())
    assert_array_almost_equal(neuron.section(0).points,
                              [[0, 0, 0], [0, 10, 0]])
    assert_array_almost_equal(neuron.section(1).points,
                              [[0, 10, 0], [-10, 10, 0]])

    neuron = Morphology(SIMPLE_PATH)
    np.random.seed(0)
    scale_morphology(neuron, ScaleParameters(mean=2, std=0.5), ScaleParameters())
    assert_array_almost_equal(neuron.section(0).points,
                              [[0, 0, 0], [0., 9.621607, 0.]])
    assert_array_almost_equal(neuron.section(1).points,
                              [[0., 9.621607, 0.], [-11.000393, 9.621607, 0.]])

    # Scaling only on Y axis
    neuron = Morphology(SIMPLE_PATH)
    scale_morphology(neuron, ScaleParameters(mean=2, axis=1), ScaleParameters())
    assert_array_almost_equal(neuron.section(0).points,
                              [[0, 0, 0], [0, 10, 0]])
    assert_array_almost_equal(neuron.section(1).points,
                              [[0, 10, 0], [-5, 10, 0]])



def test_principal_direction():
    neuron = Morphology(SIMPLE_PATH)
    assert_array_almost_equal(_principal_direction(neuron.section(0)), [0.998492, -0.0549, 0.])

    neuron.section(1).points = [[0, 0, 0], [1, 1, 0], [-1, 1, 0]]
    assert_array_equal(_principal_direction(neuron.section(1)), [1, 0, 0])

def test_create_clones():
    with TemporaryDirectory('test-create-clones') as folder:
        paths = create_clones(NEURON_PATH, folder, 1,
                              RotationParameters(30, 0, 5),
                              ScaleParameters(),
                              ScaleParameters())
        assert_equal(len(os.listdir(folder)), 1)
        assert_equal(len(paths), 1)
        out = os.listdir(folder)[0]

        actual = ImmutMorphology(os.path.join(folder, out))
        expected = ImmutMorphology(path('neuron_rotation_30_degree.swc'))
        ok_(not diff(actual, expected, atol=1e-3))

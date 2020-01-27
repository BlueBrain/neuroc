'''Module to produce clones of a morphology by jittering it'''

import attr
import numpy as np
from numpy.random import normal as normal_law
from sklearn.decomposition import PCA
from scipy.spatial.transform import Rotation

from morph_tool.transform import rotate, translate
from morphio.mut import Morphology, Section


@attr.s
class ScaleParameters:
    '''
    The scaling parameters.

    The scaling factor will be sampled upon a Normal law
    '''

    #: The mean value of the Normal law.
    #: 1 means no scaling, 2 means double the size of the object
    mean = attr.ib(type=float, default=1)

    #: The standard deviation of the Normal law.
    std = attr.ib(type=float, default=0)

    #: The axis on which to perform the scaling. If None, scaling is performed on all axes.
    axis = attr.ib(type=int, default=None)


@attr.s
class RotationParameters:
    '''
    The rotation parameters.

    Default values are from:
    https://bbpcode.epfl.ch/browse/code/platform/BlueJitterSDK/tree/apps/MorphClone.cpp#n33
    '''
    #: Mean angle (in degree)
    mean_angle = attr.ib(type=float, default=0.0)
    #: Standard deviation angle (in degree)
    std_angle = attr.ib(type=float, default=0.0)
    numberpoint = attr.ib(type=int, default=5.0)


def _principal_direction(section: Section):
    '''Get the principal direction for the cloud points made of the normed
    directions between the start of the section and all other points from this section or
    its descendent sections'''
    p0 = section.points[0]

    # Removing all first section points because they are a duplicate of the parent section last
    # point
    points = np.vstack([descendant_section.points[1:] - p0
                        for descendant_section in section.iter()])
    directions_normed = (points.T / np.linalg.norm(points, axis=1)).T
    pca = PCA(n_components=1)
    pca.fit(directions_normed)
    return pca.components_[0]


def _recursive_rotational_jitter(section: Section, piecenumber: int,
                                 angle_mean: float, angle_std: float):
    '''Rotate a section and its descendent sections

    Many rotations are applied. Each time, the rotation is applied to the section and its
    descendent sections. Similarly to how a Rubik's cube is solved:

    First applies a rotation on a leaf section
    Second applies a rotation on a leaf section and its parent section
    And so on and so on

    Args:
        section (morphio.mut.Section): the section to rotate (descendents are rotated as well)
        piecenumber (int): the number of points of the parent section to consider when making a
            rotation around the parent section
        angle_params (Dict[str, float]): the mean and std of the normal law used to sample the
            rotation angle
    '''
    for child in section.children:
        _recursive_rotational_jitter(child, piecenumber, angle_mean, angle_std)

    if section.is_root:
        return

    if not section.children:
        parent = section.parent
        direction = parent.points[-1] - parent.points[-int(min(piecenumber, len(parent.points)))]
    else:
        direction = _principal_direction(section)

    direction /= np.linalg.norm(direction)
    theta = np.random.normal(angle_mean, angle_std) * np.pi / 180.
    matrix = Rotation.from_rotvec(theta * direction).as_dcm()
    rotate(section, matrix, origin=section.points[0])


def rotational_jitter(neuron: Morphology, params: RotationParameters):
    '''Jitter sections by rotating them

    Args:
        neuron (morphio.mut.Morphology): the neuron
        params (Parameters): the parameters
    '''
    for root in neuron.root_sections:
        _recursive_rotational_jitter(root, params.numberpoint, params.mean_angle, params.std_angle)


def _segment_vectors(section: Section, prepend_null_vector: bool = False):
    '''Returns the segments of the section

    Args:
        prepend_null_vector (bool): if True, prepend [0, 0, 0] to the returned list of vectors
    '''
    vectors = np.diff(section.points, axis=0)
    if prepend_null_vector:
        return np.append([[0, 0, 0]], vectors, axis=0)
    return vectors


def _clip_minimum_scaling(scalings: float):
    '''Limit minimum scaling to 1% of original segment length'''
    return np.clip(scalings, a_min=0.01, a_max=None)


def _broadcast(scaling_factors, axis):
    '''Broadcast scaling_factor so that it applies to the specified axis.

    If axis is None, the scaling_factor is applied to all axes
    '''
    if axis is None:
        scaling_factors = np.repeat(scaling_factors[np.newaxis], 3, axis=0).T
    else:
        ones = np.ones((len(scaling_factors), 3))
        ones[:, axis] = scaling_factors.T
        scaling_factors = ones
    return scaling_factors


def _recursive_scale(section: Section,
                     segment_scaling: ScaleParameters,
                     section_scaling: ScaleParameters):
    for child in section.children:
        _recursive_scale(child, segment_scaling, section_scaling)

    # 1) Apply scaling jitter segment by segment (each segment has a different scaling factor)
    vectors = _segment_vectors(section, prepend_null_vector=True)
    scaling_factors = normal_law(segment_scaling.mean, segment_scaling.std, size=len(vectors))
    scaling_factors = _clip_minimum_scaling(scaling_factors)
    scaling_factors = _broadcast(scaling_factors, segment_scaling.axis)
    vectors = np.multiply(scaling_factors, vectors)
    cumulative_vectors = np.cumsum(vectors, axis=0)

    # 2) Apply scaling jitter at section level
    section_scaling_factor = normal_law(section_scaling.mean, section_scaling.std)
    section_scaling_factor = _clip_minimum_scaling(section_scaling_factor)
    if section_scaling.axis is None:
        cumulative_vectors *= section_scaling_factor
    else:
        cumulative_vectors[:, section_scaling.axis] *= section_scaling_factor

    new_points = section.points[0] + cumulative_vectors
    children_translation = new_points[-1] - section.points[-1]

    for child in section.children:
        translate(child, children_translation)

    section.points = new_points


def scale_morphology(neuron: Morphology,
                     segment_scaling: ScaleParameters,
                     section_scaling: ScaleParameters):
    '''
    Scale a morphology.

    The scaling is performed at 2 levels: section and segment.

    The segment scaling scales each segment of the morphology with a
    different realization of the normal law.

    The sectoin scaling scales all the segments of the same section with
    the same realization of the normal law.

    Args:
        neuron (morphio.mut.Morphology): the morphology to scale
        segment_scaling (ScaleParameters): the segment by segment specific parameters
        section_scaling (ScaleParameters): the section by section specific parameters
    '''
    for root in neuron.root_sections:
        _recursive_scale(root, segment_scaling, section_scaling)


def iter_clones(filename: str, nclones: int,
                rotation_params: RotationParameters,
                segment_scaling: ScaleParameters,
                section_scaling: ScaleParameters,
                seed: int = None):
    '''Yields 'nclones' of the input morphology

    Args:
        filename (str): the morphology to clone
        rotation_params (RotationParameters): the rotation parameters
        segment_scaling (ScaleParameters): the segment by segment specific parameters
        section_scaling (ScaleParameters): the section by section specific parameters
        seed (int): the numpy.random seed

    Yields:
        morphio.mut.Morphology clones
    '''
    if seed is not None:
        np.random.seed(seed)

    for _ in range(nclones):
        neuron = Morphology(filename)
        rotational_jitter(neuron, rotation_params)
        scale_morphology(neuron, segment_scaling, section_scaling)
        yield neuron

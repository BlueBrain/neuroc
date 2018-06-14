'''A module to perform axon shrinking '''
from __future__ import print_function

import operator
import os
import xml.etree.ElementTree

import numpy as np
from morphio import PointLevel, SectionType, set_maximum_warnings, upstream
from morphio.mut import Morphology
from tqdm import tqdm

from neuroc.utils.xml import read_placement_rules, update_rule


class NoSectionToCut(Exception):
    '''Custom exception'''


class NoAxon(Exception):
    '''Custom exception'''

    def __str__(self):
        return 'Neuron has no axon'


class TooManyAxons(Exception):
    '''Custom exception'''

    def __str__(self):
        return 'Neuron has more than one axon'


class NoAxonAnnotation(Exception):
    '''Custom exception'''

    def __str__(self):
        return 'No axon annotation'


def section_path_lengths(neuron, neurite):
    '''Return a dict {section: pathlength from the soma}'''
    dist = {s.id: np.linalg.norm(np.diff(s.points, axis=0), axis=1).sum()
            for s in neuron.iter(neurite)}

    return {section.id: sum(dist[upstream.id]
                            for upstream in neuron.iter(section, upstream))
            for section in neuron.iter(neurite)}


def cut_and_graft_coordinate(rules):
    '''Get the coordinate where to start cutting and grafting from the annotations'''
    upward = rules['dendrite']['y_min'] < rules['axon']['y_min']
    y_start_graft = rules['axon']['y_min' if upward else 'y_max']
    y_start_cut = rules['dendrite']['y_max' if upward else 'y_min']
    return upward, y_start_cut, y_start_graft


def get_start_cut_start_graft_sections(original, upward, y_start_cut, y_start_graft):
    '''Find section to cut and graft

    1) Find the axonal section with maximum pathlength
    2) Return a tuple (section_cut, section_graft) where section_cut (resp. section_graft)
    is the first section along the path to this section (starting at soma)
       whose coordinate is above y_start_cut (resp. y_start_graft)
    '''
    axons = [section for section in original.root_sections if section.type == SectionType.axon]

    if not axons:
        raise NoAxon
    if len(axons) > 1:
        raise TooManyAxons

    axon = axons[0]
    path_lengths = section_path_lengths(original, axon)
    idx_axon_end = max(path_lengths.items(), key=operator.itemgetter(1))[0]
    axon_end = original.section(idx_axon_end)

    main_branch_sections = reversed(list(original.iter(axon_end, upstream)))
    sections = {'start_cut': None, 'start_graft': None}
    for section in main_branch_sections:
        if not sections['start_cut'] and (upward == (section.points[-1, 1] > y_start_cut)):
            sections['start_cut'] = section
            continue

        if (sections['start_cut'] and not sections['start_graft'] and
                (upward == (section.points[-1, 1] > y_start_graft))):
            sections['start_graft'] = section
            break

    if not sections['start_cut']:
        raise NoSectionToCut('No section to cut from')

    if not sections['start_graft']:
        raise NoSectionToCut('No section to graft from')

    return sections['start_cut'], sections['start_graft']


def _y_interpolate(section, index1, index2, y):
    '''Return the point and diameter interpolated between points at
    'index1' and 'index2' and which is located at the given 'y' coordinate'''
    Y_INDEX = 1
    p1, p2 = section.points[[index1, index2], :]
    d1, d2 = section.diameters[[index1, index2]]
    frac = (y - p1[Y_INDEX]) / (p2[Y_INDEX] - p1[Y_INDEX])
    return (p1 + frac * (p2 - p1),
            d1 + frac * (d2 - d1))


def cut_section_at_plane_coord(section, y_plane, upward, cut_before):
    '''Cut the section so that it start (or end) at the plane coordinate
    Args:
        section (morphio.mut.Morphology): the section
        y_plane (int): the Y-coordinate at which to stop the section
        upward (bool): whether the section is oriented upward or downward
        cut_before (bool): whether to remove the points before or after the
            y_plane (according to upward/downward directionality)'''
    cut_condition = (section.points[:, 1] > y_plane) == upward
    cut_index = np.where(cut_condition)[0][0]

    if cut_index > 0:
        # Trim points before (or after) cut plane
        new_point, new_diameter = _y_interpolate(section, cut_index - 1, cut_index, y_plane)

        remaining = slice(cut_index, None) if cut_before else slice(cut_index)
        section.points = np.append([new_point] if cut_before else section.points[remaining],
                                   section.points[remaining] if cut_before else [new_point],
                                   axis=0)

        section.diameters = np.append(
            [new_diameter] if cut_before else section.diameters[remaining],
            section.diameters[remaining] if cut_before else [new_diameter],
            axis=0)


def cut_branch(neuron, upward, start_cut, y_start_cut):
    '''Delete all descendents of section start_cut and trim it so
    that it stops at y_start_cut'''

    start_cut = neuron.section(start_cut.id)
    cut_section_at_plane_coord(start_cut, y_start_cut, upward, cut_before=False)
    for child in neuron.children(start_cut):
        neuron.delete_section(child, recursive=True)


def add_vertical_segment(mut, start_cut, height):
    '''Append a vertical segment of height 'height' to section 'start_cut'
    and returns it'''
    # get the section cut in the new morphology
    # while start_cut passed as argument is the section from the original
    # morphology -> its start_cut.points[-1] has not been trimmed
    start_cut = mut.section(start_cut.id)

    vertical_segment = [start_cut.points[-1].tolist(),
                        (start_cut.points[-1] + [0, height, 0]).tolist()]

    return mut.append_section(start_cut,
                              PointLevel(vertical_segment,
                                         [start_cut.diameters[-1],
                                          start_cut.diameters[-1]]))


def translate(morphology, root, translation):
    '''Recursively translate a section and all its descendents'''
    for section in morphology.iter(root):
        section.points += translation


def graft_branch(new, original, upward, root, to_be_grafted, y_start_graft):
    '''Append section 'to_be_grafted' at the end of section 'root'
    Returns the distance along y by which section 'to_be_grafted' has been moved in the process
    '''

    cut_section_at_plane_coord(to_be_grafted, y_start_graft, upward, cut_before=True)

    translation = root.points[-1] - to_be_grafted.points[0]
    translate(original, to_be_grafted, translation)

    new.append_section(root,
                       to_be_grafted,
                       original)

    return translation[1]


def cut_and_graft(orig_filename, upward, y_start_cut, y_start_graft, height):
    '''Return the cut and grafted neuron

    1) Delete section start_cut (as well as all its descendents)
    2) Append at the point of deletion a vertical section of given height
    3) Append start_graft (and all its descendents) at the end of the vertical section
    '''

    original = Morphology(orig_filename)
    new_neuron = Morphology(orig_filename)

    start_cut, start_graft = get_start_cut_start_graft_sections(
        original, upward, y_start_cut, y_start_graft)

    cut_branch(new_neuron, upward, start_cut, y_start_cut)
    extra_section = add_vertical_segment(new_neuron, start_cut, height)
    y_min_diff = graft_branch(new_neuron, original, upward,
                              extra_section, start_graft, y_start_graft)

    return new_neuron, y_min_diff


def shrink_all_heights(orig_filename, annotation_filename,
                       output_dir,
                       heights=None,
                       n_steps=None):
    '''Perform shrinking for multiple heights of the intermediate section

    Args:
        orig_filename (str): the input morphology filename
        annotation_filename (str): the input annotation filename
        output_dir (str): the output folder

        (optional) heights (list[int]): a list of heights for the intermediate segments
        (optional) n_steps (int): the number of heights to sample if 'heights' argument
        was not passed

        One of 'heights' or 'n_steps' must be passed
        '''
    xml_tree = xml.etree.ElementTree.parse(annotation_filename)
    rules = read_placement_rules(xml_tree)
    if 'axon' not in rules:
        raise NoAxonAnnotation

    upward, y_start_cut, y_start_graft = cut_and_graft_coordinate(rules)

    if not heights:
        if not n_steps:
            raise Exception("'heights' and 'n_steps' arguments can not be None at the same time")
        y = ['y_min', 'y_max']
        height_range = abs(rules['axon'][y[not upward]] - rules['dendrite'][y[upward]])
        heights = np.linspace(0, height_range, n_steps)

    for height_added in heights:
        if not upward:
            height_added *= -1
        new_neuron, y_diff = cut_and_graft(
            orig_filename, upward, y_start_cut, y_start_graft, height_added)

        write_neuron_and_rule(new_neuron,
                              os.path.join(output_dir, '{}_height_{}'.format(
                                  os.path.basename(orig_filename).split('.')[0],
                                  height_added)),
                              xml_tree,
                              rules,
                              y_diff)


def write_neuron_and_rule(new_neuron, filename, xml_tree, rules, y_diff):
    '''Write the new neuron and new annotation to disk'''
    update_rule(xml_tree.getroot(), "axon",
                {'y_min': str(rules['axon']['y_min'] + y_diff),
                 'y_max': str(rules['axon']['y_max'] + y_diff)})

    new_neuron.write(filename + '.h5')
    xml_tree.write(filename + '.xml')


def run(file_dir, annotation_dir, output_dir, n_steps, heights):
    '''Perform the shrinking of all neurons in the file_dir'''
    set_maximum_warnings(0)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    errors = list()
    for f in tqdm(os.listdir(file_dir)):
        try:
            shrink_all_heights(os.path.join(file_dir, f),
                               os.path.join(annotation_dir, f.replace('.h5', '.xml')),
                               output_dir,
                               heights,
                               n_steps)
        except (NoSectionToCut, NoAxon, NoAxonAnnotation) as e:
            errors.append((f, str(e)))

    print('Done')
    if errors:
        print('\nThe following files could not be processed:')
        for name, error in errors:
            print('{}: {}'.format(name, error))

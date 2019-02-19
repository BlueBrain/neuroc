'''A module to perform axon shrinking '''
from __future__ import print_function

import operator
from operator import lt, gt
import os
import xml.etree.ElementTree

import numpy as np
from tqdm import tqdm
from morphio import PointLevel, SectionType, set_maximum_warnings, upstream
from morphio.mut import Morphology

from neuroc.utils.xml import read_placement_rules, update_rule

X, Y, Z = 0, 1, 2


class ShrinkError(Exception):
    '''Error related to the axon shrinker'''


def section_path_lengths(neurite):
    '''Return a dict {section: pathlength from the soma}'''
    dist = {s.id: np.linalg.norm(np.diff(s.points, axis=0), axis=1).sum()
            for s in neurite.iter()}

    return {section.id: sum(dist[upstream.id]
                            for upstream in section.iter(upstream))
            for section in neurite.iter()}


def cut_and_graft_coordinate(rules):
    '''Get the coordinate where to start cutting and grafting from the annotations'''
    upward = rules['dendrite']['y_min'] < rules['axon']['y_min']
    y_start_graft = rules['axon']['y_min' if upward else 'y_max']
    y_start_cut = rules['dendrite']['y_max' if upward else 'y_min']
    return upward, y_start_cut, y_start_graft


def get_axon(neuron):
    '''Check the neuron has exactly one axon and returns it'''
    axons = [section for section in neuron.root_sections if section.type == SectionType.axon]

    if not axons:
        raise ShrinkError('No axon')
    if len(axons) > 1:
        raise ShrinkError('Too many axons')

    return axons[0]


def get_main_branch_sections(neuron, root):
    '''Returns sections belonging to the longest path ordered from the root to the tip'''
    path_lengths = section_path_lengths(root)
    # index of the section with the longest path length
    idx_root_end = max(path_lengths.items(), key=operator.itemgetter(1))[0]
    root_end = neuron.section(idx_root_end)
    main_branch_sections = list(reversed(list(root_end.iter(upstream))))
    return main_branch_sections


def get_start_cut_start_graft_sections(original, upward, y_start_cut, y_start_graft):
    '''Find section to cut and graft

    1) Find the axonal section with maximum pathlength
    2) Return a tuple (section_cut, section_graft) where section_cut (resp. section_graft)
    is the first section along the path to this section (starting at soma)
       whose coordinate is above y_start_cut (resp. y_start_graft)
    '''
    axon = get_axon(original)
    main_branch_sections = get_main_branch_sections(original, axon)
    sections = {'start_cut': None, 'start_graft': None}
    for section in main_branch_sections:
        if not sections['start_cut']:
            if upward and (np.max(section.points[:, Y]) > y_start_cut):
                sections['start_cut'] = section
            if not upward and (np.min(section.points[:, Y]) < y_start_cut):
                sections['start_cut'] = section

        if sections['start_cut'] and not sections['start_graft']:
            if upward and (np.max(section.points[:, Y]) > y_start_graft):
                sections['start_graft'] = section
                break

            if not upward and (np.min(section.points[:, Y]) < y_start_graft):
                sections['start_graft'] = section
                break

    if not sections['start_cut']:
        raise ShrinkError('No section to cut from')

    if not sections['start_graft']:
        raise ShrinkError('No section to graft from')

    return sections['start_cut'], sections['start_graft']


def _y_interpolate(section, index1, index2, y):
    '''Return the point and diameter interpolated between points at
    'index1' and 'index2' and which is located at the given 'y' coordinate'''
    p1, p2 = section.points[[index1, index2], :]
    d1, d2 = section.diameters[[index1, index2]]
    frac = (y - p1[Y]) / (p2[Y] - p1[Y])
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
    cut_condition = (section.points[:, Y] > y_plane) == upward
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
    for child in start_cut.children:
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

    return start_cut.append_section(PointLevel(vertical_segment,
                                               [start_cut.diameters[-1],
                                                start_cut.diameters[-1]]))


def translate(root, translation):
    '''Recursively translate a section and all its descendents'''
    for section in root.iter():
        section.points += translation


def graft_branch(upward, root, to_be_grafted, y_start_graft):
    '''Append section 'to_be_grafted' at the end of section 'root'
    Returns the distance along y by which section 'to_be_grafted' has been moved in the process
    '''

    cut_section_at_plane_coord(to_be_grafted, y_start_graft, upward, cut_before=True)

    translation = root.points[-1] - to_be_grafted.points[0]
    translate(to_be_grafted, translation)

    root.append_section(to_be_grafted, recursive=True)

    return translation[1]


def cut_axon_end(neuron, y_cut):
    '''Cut axon end

    Remove points above y_cut if y_cut is above first axon point else below y_cut

    Args:
        neuron (str|morphio.mut.Morphology|morphio.Morphology):
            a morphio neuron or a neuron filename
        y_cut (float): the Y coordinate at which to cut the axon
    '''

    neuron = Morphology(neuron)
    axon = get_axon(neuron)
    upward = axon.points[0, Y] < y_cut
    main_branch_sections = get_main_branch_sections(neuron, axon)

    for section in main_branch_sections:
        op = gt if upward else lt
        idx = np.where(op(section.points, y_cut))[0]
        if idx.size:
            cut_section_at_plane_coord(section, y_cut, upward, False)
            for child in section.children:
                neuron.delete_section(child, recursive=True)
            break
    return neuron


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
    y_min_diff = graft_branch(upward, extra_section, start_graft, y_start_graft)

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
        raise ShrinkError('No axon annotation')

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
    new_neuron.write(filename + '.asc')
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
        except ShrinkError as e:
            errors.append((f, str(e)))

    print('Done')
    if errors:
        print('\nThe following files could not be processed:')
        for name, error in errors:
            print('{}: {}'.format(name, error))

'''
Module for the scaling of rat cells to human cells dimensions.
https://bbpteam.epfl.ch/project/issues/browse/IHNM-6

Scale rat cells to human cells diamensions

1) Human and rat mtypes are grouped together according to the mapping
   in MTYPE_MAPPING_FILE
2) For each group the average among all cells of the same group is computed for
   the following features:
   - standard deviation of dendritic point along Y
   - standard deviation of the radial coordinate in the XZ plane for dendritic points
   - averaged diameters of dendritic points
3) Use the ratio of the human feature to rat feature to scale rat morphologies:
   - along Y
   - along XZ
   - scale the diameters
'''

import os
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Callable, Iterable, List

import numpy as np
import yaml
from morph_tool.utils import iter_morphology_files
from morphology_repair_workflow.morphdb import MorphologyDB
from neurom import COLS, NeuriteType, iter_neurites, load_neuron
from neurom.core._neuron import Neuron
from tqdm import tqdm

from morphio.mut import Morphology
from neuroc.jitter import ScaleParameters, scale_morphology


def validate_folders(human_folder: Path, rat_folder: Path):
    '''Raises if HUMAN_FOLDER or RAT_FOLDER don't have the expected content.'''
    for path in human_folder.iterdir():
        if path.name.upper() not in {'L1', 'L2', 'L3', 'L4', 'L5', 'L6'}:
            raise ValueError('The content of the human folder: {} should only be sub-folders '
                             'with a layer name\n'
                             '(aka L1, L2, L3, etc.)\nFound: {}'.format(human_folder,
                                                                        path.name))

    neurondb = rat_folder / 'neuronDB.xml'
    if not neurondb.exists():
        raise ValueError('neuronDB.xml not found in the rat folder: {}'.format(rat_folder))


def _find_filepath(input_path: Path, stem: str):
    '''Find a morphology file with the given stem in input_path.'''
    for extension in ['.asc', '.swc', '.h5', '.ASC', '.SWC', '.H5']:
        path = Path(input_path, stem + extension)
        if path.exists():
            return path
    return None


def morphology_filenames(folder: Path):
    '''Returns a dictionary of (mtype, morphology filenames)

    from all morphologies found in the folder neuronDB.xml file
    '''
    folder = Path(folder)
    db = MorphologyDB.load_neurondb_xml(folder / 'neuronDB.xml')
    mtypes = defaultdict(list)
    for morphology in iter(db):
        path = _find_filepath(folder, morphology.name)
        if path:
            mtypes[morphology.mtype].append(path)

    return dict(mtypes)


def morphologies_in_layer(folder: Path, layer: str):
    '''Returns a list of all morphologies in the corresponding layer

    from all morphologies found in the folder neuronDB.xml file
    '''
    filename = 'neuronDB.xml'
    db = MorphologyDB.load_neurondb_xml(os.path.join(folder, filename))

    paths = list()
    for morphology in iter(db):
        path = _find_filepath(folder, morphology.name)
        if path and morphology.layer.lower() == layer[1:].lower():
            paths.append(path)

    return paths


def not_axon(neurite):
    '''Returns true if the neurite type is not an axon.'''
    return neurite.type != NeuriteType.axon


def dendritic_points(neuron: Neuron):
    '''Returns a list of all points belonging to a dendrite.'''
    return np.vstack([neurite.points[:, COLS.XYZ]
                      for neurite in iter_neurites(neuron, filt=not_axon)])


def dendritic_diameter(neuron: Neuron):
    '''Get the dendritic diameter
    '''
    radii = np.vstack([neurite.points[:, COLS.R]
                       for neurite in iter_neurites(neuron, filt=not_axon)])
    return radii.mean() * 2.


def scaling_factors(human_paths: Iterable[str],
                    rat_paths: Iterable[str],
                    funcs: List[Callable[[Neuron], float]]):
    '''Returns the list of scaling factors

    Args:
        human_paths: paths to human morphologies
        rat_paths: paths to rat morphologies
        funcs: a list of functions. For each function, a scaling factor will be computed
            by taking the ratio of the averaged value of the function applied
            to human and rat morphologies
    '''
    def means(paths):
        return np.mean([[func(morph) for func in funcs]
                        for morph in map(load_neuron, paths)],
                       axis=0)
    return means(human_paths) / means(rat_paths)


def dendritic_y_std(neuron: Neuron):
    '''Get the standard deviation of the Y coordinate for dendritic points'''
    return dendritic_points(neuron)[:, COLS.Y].std()


def dendritice_radial_std(neuron: Neuron):
    '''Get the standard deviation of the radial coordinate in the XZ plane for dendritic points'''
    radial_coord = np.linalg.norm(dendritic_points(neuron)[:, COLS.XZ],
                                  axis=1)
    return radial_coord.std()


def human_cells_per_layer_mtype(folder: Path):
    '''returns a Dict[mtype, path] of all the human cells

    FOLDER should be composed of sub-folders with the name of the layer.
    Mtype is take from the first part of the cell filename

    Example:
    - folder/L1/DAC_neuron-bla-bla.swc
    - folder/L2/AC_some-neuron.h5
    - folder/L3/BTF_another-neuron.ASC
    '''
    cells = dict()
    for layer_folder in folder.iterdir():
        layer = layer_folder.stem
        cells[layer] = defaultdict(list)
        for cell in iter_morphology_files(layer_folder):
            mtype = cell.stem.split('_')[0]
            cells[layer][mtype].append(cell)
    return cells


def mtype_name(layer: str, mtype: str):
    '''Returns the standardized name of the mtype for a given layer.'''
    L23_mtypes = {'LBC', 'BP', 'BTC', 'CHC', 'DBC', 'LBC', 'MC', 'NBC', 'NGC', 'SBC'}
    if layer.upper() in {'L2', 'L3'} and mtype in L23_mtypes:
        return 'L23_' + mtype
    return layer + '_' + mtype


def mtype_matcher(human_folder: Path,
                  rat_folder: Path,
                  mtype_mapping: Path):
    '''Group human and rat cells by equivalence of mtype.

    Human and rat cells have different mtypes but can be grouped together by using
    a mtype mapping.

    Returns a list of tuples (human cells, rat cells) whose mtypes can be considered
    as equivalent.
    '''
    zipped = list()

    with mtype_mapping.open() as file_:
        mtype_mapping = yaml.load(file_, Loader=yaml.FullLoader)

    rat_morphologies = morphology_filenames(rat_folder)
    human_cells_dict = human_cells_per_layer_mtype(human_folder)

    for layer in mtype_mapping:
        for human_mtype, rat_mtypes in mtype_mapping[layer].items():
            rat_cells = list()

            if human_mtype.lower() == 'all':
                human_cells = list(iter_morphology_files(Path(human_folder, layer)))
            else:
                human_cells = human_cells_dict[layer][human_mtype]

            if len(rat_mtypes) == 1 and rat_mtypes[0].lower() == 'all':
                zipped.append((human_cells, morphologies_in_layer(rat_folder, layer)))
                continue

            for mtype in rat_mtypes:
                if mtype.lower() == 'all':
                    raise ValueError('The human -> rat mtype mapping is: {} -> {}\n'
                                     'If you specify "all", it should be the only element'
                                     ' of the mapping'.format(human_mtype, rat_mtypes))

                if mtype_name(layer, mtype) not in rat_morphologies:
                    warnings.warn('{} not found among available rat cells'.format(
                        mtype_name(layer, mtype)))
                    continue
                rat_cells.extend(rat_morphologies[mtype_name(layer, mtype)])

            if human_cells or rat_cells:
                zipped.append((human_cells, rat_cells))

    return zipped


def iter_scaling_and_rat(human_folder: Path,
                         rat_folder: Path,
                         mtype_mapping_file: Path,
                         funcs: List[Callable[[Neuron], float]]):
    '''Yiels rat cells and their scaling factors.

    Args:
        human_folder: path to human folder
        rat_folder: path to rat folder
        mtype_mapping_file: path to YAML file containing the mtype mapping
        funcs: list of functions returning a neuron feature. The average of the
            feature among groups (human/rat) of mapped cells will be used to compute
            the scaling factor
    '''
    for humans, rats in mtype_matcher(human_folder, rat_folder, mtype_mapping_file):
        if not humans:
            continue
        if not rats:
            warnings.warn('No matching rat morphologies for human mtypes: {}\n'
                          'Skipping the scaling'.format(humans))
            continue

        for rat in rats:
            yield rat, scaling_factors(humans, rats, funcs)


def scale_diameter(neuron: Morphology,
                   scaling_factor: float):
    '''Scale a neuron by a constant scaling factor:

    Args:
        neuron: a neuron
        scaling_factor: the diameter scaling factor (ie. 2 means it doubles the diameters)
    '''
    for section in neuron.iter():
        section.diameters = scaling_factor * section.diameters


def scale_rat_cells(human_folder: Path,
                    rat_folder: Path,
                    mtype_mapping_file: Path,
                    output_folder: Path):
    '''Scale rat cells to human cells diamensions

    HUMAN_DIR should be a dir with the following structure:
        - Must be **only** composed of sub-folders whose filename is a layer name
        - Each sub folder should be composed of morphology files whose first part of the filename
      before the '_' is considered as the **mtype**

    RAT_DIR should be a directory containing rat morphology files **and a neuronDB.xml file.

    MTYPE_MAPPING_FILE is a YAML file containing a dictionary where:
       - a key is a human mtype or **all**
       - the value is a list of rat mtypes to associate with the key. Or a list of one 'all' element


    Algorithm:
    1) Human and rat mtypes are grouped together according to the mapping
       in MTYPE_MAPPING_FILE
    2) For each group the average among all cells of the same group is computed for
       the following features:
       - standard deviation of dendritic point along Y
       - standard deviation of the radial coordinate in the XZ plane for dendritic points
       - averaged diameters of dendritic points
    3) Use the ratio of the human feature to rat feature to scale rat morphologies:
       - Use 1st feature to scale along Y
       - Use 2nd feature to scale along XZ
       - Use 3rd feature to scale the diameters

    See issue:
    https://bbpteam.epfl.ch/project/issues/browse/IHNM-6
    '''
    validate_folders(human_folder, rat_folder)
    scaling_variables = (dendritic_y_std, dendritice_radial_std, dendritic_diameter)
    rats_and_factors = iter_scaling_and_rat(human_folder,
                                            rat_folder,
                                            mtype_mapping_file,
                                            scaling_variables)
    for rat, (y_scale, xz_scale, diam_scale) in tqdm(list(rats_and_factors)):
        path = Path(rat)
        neuron = Morphology(path)
        scale_morphology(neuron,
                         ScaleParameters(),  # no segment level scaling is performed
                         ScaleParameters(mean=y_scale, axis=COLS.Y))
        scale_morphology(neuron,
                         ScaleParameters(),  # no segment level scaling is performed
                         ScaleParameters(mean=xz_scale, axis=COLS.XZ))
        scale_diameter(neuron, diam_scale)
        name_extension = '_-_Y-Scale_{}_-_XZ-Scale_{}_-_Diam-Scale_{}.h5'.format(
            y_scale, xz_scale, diam_scale)
        neuron.write(Path(output_folder, path.stem + name_extension))

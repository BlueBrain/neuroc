# SPDX-License-Identifier: Apache-2.0
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

import logging
from pathlib import Path
from typing import Callable, Iterable, List, Tuple

import numpy as np
import pandas as pd
import yaml
from neurom import COLS, NeuriteType, iter_neurites, load_morphology
from neurom.core import Morphology
from tqdm import tqdm

L = logging.getLogger('neuroc')


def neurondb_dataframe(filename: Path) -> pd.DataFrame:
    '''Returns a DataFrame: [name, layer, mtype]

    Args:
        filename: the neurondb.dat file
    '''
    df = pd.read_csv(filename, sep=' ', names=['name', 'layer', 'mtype'], index_col=False)
    df.layer = df.layer.astype('str')
    return df


def _find_filepath(input_path: Path, stem: str):
    '''Find a morphology file with the given stem in input_path.'''
    children = set(input_path.iterdir())
    for extension in ['.asc', '.swc', '.h5', '.ASC', '.SWC', '.H5']:
        path = Path(input_path, stem + extension)
        # check in children for case sensitivity
        if path in children:
            return path.resolve()
    raise ValueError(f'{stem} not found in {input_path}')


def not_axon(neurite):
    '''Returns true if the neurite type is not an axon.'''
    return neurite.type != NeuriteType.axon


def dendritic_points(neuron: Morphology):
    '''Returns a list of all points belonging to a dendrite.'''
    return np.vstack([neurite.points[:, COLS.XYZ]
                      for neurite in iter_neurites(neuron, filt=not_axon)])


def dendritic_diameter(neuron: Morphology):
    '''Get the dendritic diameter
    '''
    radii = np.hstack([neurite.points[:, COLS.R]
                       for neurite in iter_neurites(neuron, filt=not_axon)])
    return radii.mean() * 2.


def scaling_factors(human_paths: Iterable[str],
                    rat_paths: Iterable[str],
                    funcs: Iterable[Callable[[Morphology], float]]):
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
                        for morph in map(load_morphology, paths)],
                       axis=0)
    return means(human_paths) / means(rat_paths)


def dendritic_y_std(neuron: Morphology):
    '''Get the standard deviation of the Y coordinate for dendritic points'''
    return dendritic_points(neuron)[:, COLS.Y].std()


def dendritice_radial_std(neuron: Morphology):
    '''Get the standard deviation of the radial coordinate in the XZ plane for dendritic points'''
    radial_coord = np.linalg.norm(dendritic_points(neuron)[:, COLS.XZ],
                                  axis=1)
    return radial_coord.std()


def extended_neurondb(morphology_neurondb: Path) -> pd.DataFrame:
    '''Returns a DataFrame [layer, mtype, name, path]'''
    df = neurondb_dataframe(morphology_neurondb)
    folder = morphology_neurondb.parent.joinpath('h5')
    df['path'] = [_find_filepath(folder, name) for name in df.name]
    return df


def mtype_matcher(human_neurondb: Path,
                  rat_neurondb: Path,
                  mtype_mapping: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    '''Group human and rat cells by equivalence of mtype.

    Human and rat cells have different mtypes but can be grouped together by using
    a mtype mapping.

    Yields 2-tuples of dataframes of (human cells, rat cells)
    whose mtypes can be considered as equivalent.
    '''
    with mtype_mapping.open() as file_:
        mtype_mapping = yaml.load(file_, Loader=yaml.FullLoader)

    df_rat = extended_neurondb(rat_neurondb)
    df_human = extended_neurondb(human_neurondb)

    for human_mtype, rat_mtypes in mtype_mapping.items():
        rats = df_rat[df_rat.mtype.isin(rat_mtypes)]
        assert not rats.empty, (f'The following {human_mtype} (human) -> {rat_mtypes} (rat)'
                                ' did not return any rat cells')
        humans = df_human[df_human.mtype == human_mtype]
        assert not humans.empty, (f'No morphology found for human mtype: {human_mtype}')
        yield humans, rats


def iter_scaling_and_rat(
    human_neurondb: Path,
    rat_neurondb: Path,
    mtype_mapping_file: Path,
    funcs: Iterable[Callable[[Morphology], float]]
) -> Tuple[str, str, str, Path, List[float]]:
    '''Yields a tuple (human layer, rat mtype, rat layer, rat path, scaling factors)

    Args:
        human_neurondb: path to human neurondb
        rat_neurondb: path to rat neurondb
        mtype_mapping_file: path to YAML file containing the mtype mapping
        funcs: list of functions returning a neuron feature. The average of the
            feature among groups (human/rat) of mapped cells will be used to compute
            the scaling factor
    '''
    sequence = list(mtype_matcher(human_neurondb, rat_neurondb, mtype_mapping_file))
    for df_human, df_rat in tqdm(sequence):
        factors = scaling_factors(df_human.path, df_rat.path, funcs)

        human_mtype = df_human.mtype.iloc[0]
        for rat_mtype, rat_layer, rat_path in zip(df_rat.mtype, df_rat.layer, df_rat.path):
            yield human_mtype, rat_mtype, rat_layer, rat_path, factors


def scale_diameter(neuron: Morphology,
                   scaling_factor: float) -> None:
    '''In-place scaling of a neuron diameters by a constant scaling factor.

    Args:
        neuron: a neuron
        scaling_factor: the diameter scaling factor (ie. 2 means it doubles the diameters)
    '''
    for section in neuron.iter():
        section.diameters = scaling_factor * section.diameters


def scale_coordinates(neuron: Morphology, scaling_factor: float, axis: COLS) -> None:
    '''In-place scaling of a neuron neurites coordinate by a constant scaling factor.

    Args:
        neuron: a neuron
        scaling_factor: the diameter scaling factor (ie. 2 means it doubles the coordinates)
        axis: the NeuroM axis (or axes) on which to perform the scaling
    '''
    for section in neuron.iter():
        points = np.copy(section.points)
        points[:, axis] *= scaling_factor
        section.points = points


def scale_one_cell(path: Path, y_scale: float, xz_scale: float, diam_scale: float) -> Morphology:
    '''Scale the coordinates and diameters of a cell.

    Args:
        path: a neuron path
        y_scale: Y scaling factor
        xz_scale: XZ scaling factor
        diam_scale: diameter scaling factor

    Returns:
        a scaled morphology
    '''
    neuron = Morphology(path)
    scale_coordinates(neuron, y_scale, COLS.Y)
    scale_coordinates(neuron, xz_scale, COLS.XZ)
    scale_diameter(neuron, diam_scale)
    return neuron


def scale_all_cells(human_neurondb: Path,
                    rat_neurondb: Path,
                    mtype_mapping_file: Path,
                    output_folder: Path) -> None:
    '''Scale rat cells to human cells diamensions

    Args:
        human_neurondb: the human neurondb filename
        rat_neurondb: the rat neurondb filename
        mtype_mapping_file: the YAML mapping HUMAN mtype to RAT mtypes
            It must be a dictionary
                key: human mtype
                value: list of corresponding rat mtypes
        output_folder: the output folder

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
    scaling_functions = (dendritic_y_std, dendritice_radial_std, dendritic_diameter)

    L.info('Grouping human and rat cells by equivalent mtypes. This may take a while...')
    iterable = list(iter_scaling_and_rat(human_neurondb,
                                         rat_neurondb,
                                         mtype_mapping_file,
                                         scaling_functions))

    L.info('Scaling rat cells. This may take a while...')

    metadata = []

    for human_mtype, rat_mtype, rat_layer, rat_path, (y_scale, xz_scale, diam_scale) in tqdm(
            iterable):
        name = (f'{rat_path.stem}_-_Y-Scale_{y_scale}_-_XZ-Scale_{xz_scale}_-_Diam-'
                f'Scale_{diam_scale}.h5')
        scale_one_cell(rat_path, y_scale, xz_scale, diam_scale).write(output_folder / name)
        metadata.append([Path(name).stem, human_mtype, rat_mtype, rat_layer, rat_path.stem])

    pd.DataFrame(
        data=metadata,
        columns=['human_name', 'human_mtype', 'rat_mtype', 'rat_layer', 'original_rat']
    ).to_csv(output_folder / 'neurondb.csv', index=False)

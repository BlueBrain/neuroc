'''CLI endpoint'''
# pylint: disable=import-outside-toplevel
from pathlib import Path
import click

import morphio
from morph_tool.utils import iter_morphology_files

from neuroc.axon_shrinker.shrink import run
from neuroc.axon_shrinker.viewer import app, set_output_folder, set_input_folder


REQUIRED_DIR = click.Path(exists=True, readable=True, dir_okay=True, resolve_path=True)


@click.group()
def cli():
    '''The CLI object'''


@cli.group()
def scale():
    '''Scaling utilities'''


@scale.group()
def simple():
    '''Scale morphologies with a constant scaling factor.

    Note: it does not scale the diameter
    '''

# pylint: disable=function-redefined
@simple.command(short_help='Scale one morphology')
@click.argument('input_file', type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('output_file')
@click.option('--scaling', type=float, required=True,
              help='The scaling value')
def file(input_file, output_file, scaling):
    '''Scale a morphology with a constant scaling factor.

    Note: it does not scale the diameter
    '''
    from neuroc.jitter import ScaleParameters, scale_morphology
    neuron = morphio.mut.Morphology(input_file)
    scale_morphology(neuron,
                     segment_scaling=ScaleParameters(),
                     section_scaling=ScaleParameters(mean=scaling))
    neuron.write(output_file)


# pylint: disable=function-redefined
@simple.command(short_help='Scale all morphologies in a folder')
@click.argument('input_dir')
@click.argument('output_dir', type=click.Path(exists=True, file_okay=False, writable=True))
@click.option('--scaling', type=float, required=True,
              help='The scaling value')
def folder(input_dir, output_dir, scaling):
    '''Scale all morphologies in the folder with a constant scaling factor.

    Note: it does not scale the diameter
    '''
    from neuroc.jitter import ScaleParameters, scale_morphology
    for path in iter_morphology_files(input_dir):
        neuron = morphio.mut.Morphology(path)
        scale_morphology(neuron,
                         segment_scaling=ScaleParameters(),
                         section_scaling=ScaleParameters(mean=scaling))
        neuron.write(str(Path(output_dir, Path(path).name)))


@cli.command(short_help='Shrink an axon using a web app')
@click.argument('output_folder', type=click.Path(exists=True, file_okay=False, writable=True))
@click.argument('input_folder', type=click.Path(exists=True, file_okay=False))
def axon_shrinker_viewer(input_folder, output_folder):
    '''Open the webapp to shrink an axon manually'''
    set_input_folder(input_folder)
    set_output_folder(output_folder)
    app.run_server(debug=True)


@cli.command(short_help='Clone morphologies and splice their axon')
@click.argument('files_folder', type=click.Path(exists=True, file_okay=False))
@click.argument('annotations_folder', type=click.Path(exists=True, file_okay=False))
@click.argument('output_folder', type=click.Path(exists=True, file_okay=False, writable=True))
@click.option('--nsamples', default=10)
@click.option('--heights', default=None, multiple=True, type=int)
def axon_shrinker(files_folder, annotations_folder, output_folder, nsamples, heights):
    '''For each morphology of the FILES_FOLDER, remove the axon splice described by the
    corresponding annotation (ie. located between the end of the dendritic annotation
    and the start of the axonal annotation) and replace it by an intermediate vertical segment.

    For each input morphology, the length of the replaced segment is either determined by
    values in the list argument HEIGHTS if provided, else heights will be sampled from 0 to
    the length of initially spliced segment. In this case, the number of samples can be passed
    with the NSAMPLES argument (default to 10).

    example:

    \b
    neuroc axon_shrinker files_dir annotations_dir output_dir
    '''
    run(files_folder, annotations_folder, output_folder, nsamples, heights)


# pylint: disable=function-redefined
@scale.command(short_help='Scale rat cell to human cell dimensions')
@click.argument('human_dir', type=REQUIRED_DIR)
@click.argument('rat_dir', type=REQUIRED_DIR)
@click.argument('mtype_mapping', type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('output_dir', type=click.Path(exists=True, file_okay=False, writable=True))
def rat_to_human(human_dir, rat_dir, mtype_mapping, output_dir):
    '''Scale rat cells to human cells diamensions

    HUMAN_DIR should be a dir with the following structure:
        - Must be **only** composed of sub-folders whose filename is a layer name
        - Each sub folder should be composed of morphology files whose first part of the filename
      before the '_' is considered as the **mtype**

    RAT_DIR should be a directory containing rat morphology files **and a neuronDB.xml file.

    MTYPE_MAPPING_FILE is a YAML file containing a dictionary where:
       - a key is a human mtype or **all**
       - the value is a list of rat mtypes to associate with the key. Or a list of one 'all' element

    \b
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
    from neuroc.rat_to_human import scale_all_cells
    scale_all_cells(Path(human_dir),
                    Path(rat_dir),
                    Path(mtype_mapping),
                    Path(output_dir))

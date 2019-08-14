'''CLI endpoint'''
import click

import numpy as np

from neuroc.axon_shrinker.shrink import run
from neuroc.axon_shrinker.viewer import app, set_output_folder, set_input_folder

from neuroc.jitter import ScaleParameters, RotationParameters, create_clones


@click.group()
def cli():
    '''The CLI object'''


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


# pylint: disable=too-many-arguments
@cli.command(short_help='Create clones of the input morphology and jitter the neurites')
@click.argument('morphology', type=str)
@click.option('-o', '--output_folder', type=click.Path(exists=True, file_okay=False, writable=True),
              default='.', help='The output directory')
@click.option('-n', '--n-clones', help='Number of clones to generate', default=10)
@click.option('-R', '--mean-angle', help='Mean of the normal distribution of angle (in degrees)',
              default=0.)
@click.option('-r', '--std-angle',
              help='Standard Deviation of the normal distribution of angle'
              ' (in degrees)', default=0.)
@click.option('-S', '--mean-segment-scaling', help='Mean of the normal distribution for the scaling'
              ' factor segment-wise (in %, ie. 100 = times 2 scaling)', default=0.)
@click.option('-s', '--std-segment-scaling',
              help='Standard deviation of the normal distribution for the scaling factor'
              ' segment-wise (in %)', default=0.)
@click.option('-B', '--mean-section-scaling', help='Mean of the normal distribution for the scaling'
              ' factor section-wise (in %)', default=0.)
@click.option('-b', '--std-section-scaling', help='Standard deviation of the normal distribution'
              ' for the scaling factor section-wise (in %)', default=0.)
@click.option('--numberpoint', default=5,
              help='Number of points used to define the rotation direction when ROTALGO="pat"')
@click.option('-d', '--seed', default=0, help='The numpy random seed')
def clone_army(morphology, output_folder, n_clones, mean_angle, std_angle, mean_segment_scaling,
               std_segment_scaling, mean_section_scaling, std_section_scaling,
               numberpoint, seed):
    """Create clones of the input morphology

    Some jitter is added to make them all different. There are two kinds of jitters:

    - rotation: each section around its parent axis or around the PCA of all descendant points
    - scaling: individual segment lengths as well as section lengths are jittered

    Example:

    \b
    neuroc clone-army input-morph.asc output-dir --rmean 30 --rstd 10
    """
    np.random.seed(seed)
    segment_scaling = ScaleParameters(mean=mean_segment_scaling / 100.,
                                      std=std_segment_scaling / 100.)
    section_scaling = ScaleParameters(mean=mean_section_scaling / 100.,
                                      std=std_section_scaling / 100.)
    rotation_params = RotationParameters(mean_angle=mean_angle,
                                         std_angle=std_angle,
                                         numberpoint=numberpoint)

    create_clones(morphology, output_folder, n_clones, rotation_params,
                  segment_scaling, section_scaling)

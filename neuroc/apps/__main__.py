import click

from neuroc.axon_shrinker.shrink import run


@click.group()
def cli():
    pass


@cli.command(short_help='Clone morphologies and splice their axon')
@click.argument('files_folder')
@click.argument('annotations_folder')
@click.argument('output_folder')
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

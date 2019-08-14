from tempfile import TemporaryDirectory
from pathlib import Path
from os.path import dirname, join as joinp
from nose.tools import assert_equal
from click.testing import CliRunner
from neuroc.cli import cli


DATA_PATH = Path(Path(__file__).parent, 'data')

def path(file):
    return str(Path(DATA_PATH, file))

def test_clone_army():
    runner = CliRunner()
    with TemporaryDirectory() as output_folder:
        result = runner.invoke(cli, ['clone-army', path('simple.asc'),
                                     '-o', output_folder])
    assert_equal(result.exit_code, 0)

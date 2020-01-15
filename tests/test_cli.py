from pathlib import Path

from tempfile import TemporaryDirectory
from click.testing import CliRunner
from nose.tools import assert_equal

from neuroc.cli import cli

PATH = Path(__file__).resolve().parent / 'data'


def test_cli():
    runner = CliRunner()

    with TemporaryDirectory(prefix='test-scale-file') as folder:
        folder = Path(folder)
        result = runner.invoke(cli, ['scale', 'file',
                                     '--scaling', '2.0',
                                     str(PATH / 'simple.asc'),
                                     str(folder / 'simple-scaled.asc')])
        assert_equal(result.exit_code, 0, result.exception)
        assert_equal(len(list(folder.glob('*'))), 1)

    with TemporaryDirectory(prefix='test-scale-folder') as folder:
        folder = Path(folder)
        result = runner.invoke(cli, ['scale', 'folder',
                                     '--scaling', '2.0',
                                     str(PATH),
                                     str(folder)])
        assert_equal(result.exit_code, 0, result.exception)
        assert_equal(len(list(folder.glob('*'))),
					 len(list(PATH.glob('*'))))

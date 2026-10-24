import os
import sys

import pytest
import yaml

pytest.importorskip('plotly')
pytest.importorskip('IPython')

from conftest import REPO

sys.path.insert(0, os.path.join(REPO, 'examples'))
from be3d_local_helper import EXAMPLE_CONFIGS, localize_config, notebook_mode, prepare_run_folders  # noqa: E402


@pytest.fixture
def run_dir(tmp_path):
    run = str(tmp_path / 'BE3D_run')
    prepare_run_folders(REPO, run)
    return run


def write_yaml(path, cfg):
    with open(path, 'w') as f:
        yaml.safe_dump(cfg, f)
    return str(path)


def test_prepare_run_folders_copies_examples_without_overwriting(tmp_path):
    run = str(tmp_path / 'BE3D_run')
    dirs = prepare_run_folders(REPO, run)
    assert sorted(os.listdir(dirs['structures'])) == ['8voj_KBTBD4_HDAC1_chain_ABC.pdb', 'MORC2.pdb']
    assert 'Human_InVitro_HS294T_Apobec_D24_Input.tsv' in os.listdir(dirs['screens'])
    assert os.path.isdir(dirs['configs']) and os.path.isdir(dirs['outputs'])
    mine = os.path.join(dirs['structures'], 'MORC2.pdb')
    with open(mine, 'w') as f:
        f.write('my own structure\n')
    prepare_run_folders(REPO, run)
    assert open(mine).read() == 'my own structure\n'


@pytest.mark.parametrize('mode', sorted(EXAMPLE_CONFIGS))
def test_bundled_configs_find_all_their_inputs(run_dir, tmp_path, mode):
    src = os.path.join(REPO, 'examples', 'yaml', EXAMPLE_CONFIGS[mode])
    cfg, missing, _ = localize_config(src, str(tmp_path / 'out.yaml'), run_dir, '/repo', '/content/muscle')
    assert missing == []
    assert notebook_mode(cfg) == mode
    assert cfg['output_dir'].startswith(os.path.join(run_dir, 'outputs'))
    assert cfg['beclust3d_path'] == '/repo' and cfg['muscle_path'] == '/content/muscle'


def test_foreign_paths_are_resolved_by_file_name(run_dir, tmp_path):
    screens = os.path.join(run_dir, 'inputs', 'screens')
    open(os.path.join(screens, 'mine.tsv'), 'w').write('x\n')
    src = write_yaml(tmp_path / 'my.yaml', {
        'mode': 'monomer', 'screen_dir': '/home/someone/data', 'screens': 'mine.tsv',
        'user_pdb': '/home/someone/pdbs/MORC2.pdb', 'output_dir': '/home/someone/results/MyRun',
        'conservation': {'run': True, 'alignment_filename': '/home/someone/aln/my.align'},
    })
    # the real checkout, so the bundled examples can be told apart from the user's own files
    cfg, missing, from_examples = localize_config(src, str(tmp_path / 'my_colab.yaml'), run_dir, REPO, '/m')
    assert cfg['screen_dir'] == screens
    assert cfg['user_pdb'] == os.path.join(run_dir, 'inputs', 'structures', 'MORC2.pdb')
    assert cfg['output_dir'] == os.path.join(run_dir, 'outputs', 'MyRun')
    assert missing == [('conservation.alignment_filename', 'my.align')]
    # MORC2.pdb was found only as the bundled example of that name
    assert from_examples == [('user_pdb', cfg['user_pdb'])]
    assert notebook_mode(cfg) == 'cross_species'
    assert '/home/someone/data' in open(src).read()  # the user's own file is left alone


def test_missing_screens_are_listed_and_partners_are_handled(run_dir, tmp_path):
    src = write_yaml(tmp_path / 'blind.yaml', {
        'mode': 'blind_target', 'output_dir': '/x/out',
        'partners': [{'screen_dir': '/elsewhere', 'screens': 'YeoKBTBD4-HDAC12025-ABE-Screen-PPI.tsv, absent.tsv',
                      'muscle_path': '/usr/bin/muscle'}],
    })
    cfg, missing, _ = localize_config(src, str(tmp_path / 'blind_colab.yaml'), run_dir, '/repo', '/m')
    assert missing == [('partners[0].screens', 'absent.tsv')]
    assert cfg['partners'][0]['muscle_path'] == '/m'


def test_notebook_mode_mapping():
    assert notebook_mode({}) == 'monomer'
    assert notebook_mode({'mode': 'monomer', 'conservation': {'run': True}}) == 'cross_species'
    assert notebook_mode({'mode': 'ppi_diff'}) == 'ppi'
    assert notebook_mode({'mode': 'blind_target'}) == 'blind_target'
    with pytest.raises(ValueError, match='ppi_diff'):
        notebook_mode({'mode': 'complex'})

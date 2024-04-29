
from shutil import copy
import yaml

files = {
    #snakemake.config['configfile']: "config.yaml",
    "config.yaml": "config.yaml",
    "Snakefile": "Snakefile",
    "scripts/solve_network.py": "solve_network.py",
    "scripts/prepare_sector_network.py": "prepare_sector_network.py",
    "../pypsa-eur/config.yaml": "config.pypsaeur.yaml"
}

if __name__ == '__main__':
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake('copy_config')

    basepath = snakemake.config['summary_dir'] + '/' + snakemake.config['run'] + '/configs/'
    #print(snakemake.config['configfile'])

    for f, name in files.items():
        print(f)
        print(basepath + name)
        copy(f, basepath + name)


    with open(basepath + 'config.snakemake.yaml', 'w') as yaml_file:
        yaml.dump(
            snakemake.config,
            yaml_file,
            default_flow_style=False,
            allow_unicode=True,
            sort_keys=False
        )

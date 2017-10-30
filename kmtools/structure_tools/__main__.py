import logging

import fire

from kmtools import structure_tools


def main():
    fire.Fire({
        'interaction_dataset':
        structure_tools.interaction_dataset.generate_interaction_dataset
    })
    return 0


if __name__ == '__main__':
    import sys
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    logging.getLogger('kmbio.PDB.Atom').setLevel(logging.WARNING)
    sys.exit(main())

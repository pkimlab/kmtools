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
    sys.exit(main())

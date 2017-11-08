"""
Process ``components.cif.gz``, downloaded from https://www.wwpdb.org/data/ccd#mmcifFormat,
into a dictionary mapping modified residues to canonical residues.
"""
import fileinput
import json
import pickle
import shlex


def process_chemical_components_dictionary():
    component_dict = {}
    # Instance variables
    key = None
    value = None
    component_id = None
    component_data = {}
    value_flag = None
    for i, line in enumerate(fileinput.input()):
        if line.startswith('data_'):
            if component_id and component_data:
                component_dict[component_id] = component_data
            component_id = line[5:].strip()
            component_data = {}
        elif value_flag == 1:
            assert line.startswith(';')
            value += shlex.split(line.strip().strip(';'))[0]
            value_flag = 2
        elif value_flag == 2:
            if line.startswith(';'):
                assert not line.strip().strip(';')
                component_data[key] = value if value != '?' else None
                value_flag = None
            else:
                value += shlex.split(line.strip())[0]
        else:
            for key in ['_chem_comp.id', '_chem_comp.type', '_chem_comp.mon_nstd_parent_comp_id']:
                if line.startswith(key):
                    if len(line) > 50:
                        value = shlex.split(line.strip())[-1]
                        component_data[key] = value if value != '?' else None
                    else:
                        value = ""
                        value_flag = 1
                        break
    return component_dict


def gen_modified_to_canonical_mapping(chemical_components):
    """Generate a dictionary mapping modified residues to canonical residues."""
    modified_to_canonical = {
        key: data['_chem_comp.mon_nstd_parent_comp_id'].upper().split(',')[0]
        for key, data in chemical_components.items()
        if data['_chem_comp.mon_nstd_parent_comp_id'] is not None and
        'peptide' in data['_chem_comp.type'].lower()
    }
    return modified_to_canonical


if __name__ == '__main__':
    chemical_components = process_chemical_components_dictionary()
    with open('chemical_components.pickle', 'wb') as ofh:
        pickle.dump(chemical_components, ofh, pickle.HIGHEST_PROTOCOL)

    modified_to_canonical = gen_modified_to_canonical_mapping(chemical_components)
    with open('../kmtools/structure_tools/data/residue_mapping_to_canonical.json', 'wt') as fout:
        json.dump(modified_to_canonical, fout, indent=4, sort_keys=True)

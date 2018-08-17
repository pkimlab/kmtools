"""
Process ``components.cif.gz``, downloaded from http://www.wwpdb.org/data/ccd (exact link:
ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif) into a dictionary mapping
modified residues to canonical residues.
"""
import fileinput
import json
import pickle
import shlex

RNA = ["A", "T", "C", "G", "U"]
DNA = ["DA", "DT", "DC", "DG", "DU", "DI"]
AAA = [
    "ALA",
    "ASN",
    "VAL",
    "PRO",
    "CYS",
    "DCY",  #
    "ASP",
    "TYR",
    "LYS",
    "PHE",
    "GLY",
    "TRP",
    "LEU",
    "SER",
    "ARG",
    "THR",
    "GLU",
    "MET",
    "HIS",
    "SEC",  #
    "ILE",
    "GLN",
    "FUC",  #
    "DSN",  #
    "DTH",  #
    "PYL",  #
]


def process_chemical_components_dictionary():
    component_dict = {}
    # Instance variables
    key = None
    value = None
    component_id = None
    component_data = {}
    value_flag = None
    for i, line in enumerate(fileinput.input()):
        if line.startswith("data_"):
            if component_id and component_data:
                component_dict[component_id] = component_data
            component_id = line[5:].strip()
            component_data = {}
        elif value_flag == 1:
            assert line.startswith(";")
            value += shlex.split(line.strip().strip(";"))[0]
            value_flag = 2
        elif value_flag == 2:
            if line.startswith(";"):
                assert not line.strip().strip(";")
                component_data[key] = value if value != "?" else None
                value_flag = None
            else:
                value += shlex.split(line.strip())[0]
        else:
            for key in ["_chem_comp.id", "_chem_comp.type", "_chem_comp.mon_nstd_parent_comp_id"]:
                if line.startswith(key):
                    if len(line) > 50:
                        value = shlex.split(line.strip())[-1]
                        component_data[key] = value if value != "?" else None
                    else:
                        value = ""
                        value_flag = 1
                        break
    return component_dict


def gen_modified_to_canonical_mapping(chemical_components):
    """Generate a dictionary mapping modified residues to canonical residues."""
    modified_to_canonical = {
        key: data["_chem_comp.mon_nstd_parent_comp_id"].upper().split(",")[0]
        for key, data in chemical_components.items()
        if data["_chem_comp.mon_nstd_parent_comp_id"] is not None
        # and "linking" in data["_chem_comp.type"].lower()
    }
    return modified_to_canonical


if __name__ == "__main__":
    chemical_components = process_chemical_components_dictionary()
    with open("components.pickle", "wb") as ofh:
        pickle.dump(chemical_components, ofh, pickle.HIGHEST_PROTOCOL)

    modified_to_canonical = gen_modified_to_canonical_mapping(chemical_components)

    with open("../kmtools/structure_tools/data/rna_mapping_to_canonical.json", "wt") as fout:
        rna_mapping_to_canonical = {k: v for k, v in modified_to_canonical.items() if len(v) == 1}
        assert {v for v in rna_mapping_to_canonical.values()} == set(RNA), {
            v for v in rna_mapping_to_canonical.values()
        }
        json.dump(rna_mapping_to_canonical, fout, indent=4, sort_keys=True)

    with open("../kmtools/structure_tools/data/dna_mapping_to_canonical.json", "wt") as fout:
        dna_mapping_to_canonical = {k: v for k, v in modified_to_canonical.items() if len(v) == 2}
        assert {v for v in dna_mapping_to_canonical.values()} == set(DNA), {
            v for v in dna_mapping_to_canonical.values()
        }
        json.dump(dna_mapping_to_canonical, fout, indent=4, sort_keys=True)

    with open("../kmtools/structure_tools/data/aaa_mapping_to_canonical.json", "wt") as fout:
        aaa_mapping_to_canonical = {k: v for k, v in modified_to_canonical.items() if len(v) == 3}
        assert {v for v in aaa_mapping_to_canonical.values()} == set(AAA), {
            v for v in aaa_mapping_to_canonical.values()
        }
        json.dump(aaa_mapping_to_canonical, fout, indent=4, sort_keys=True)

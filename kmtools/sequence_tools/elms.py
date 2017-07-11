import os
import re
import pandas as pd


# To improve. Motif definitios repositories
absolute_path = os.path.dirname(__file__)


# ELM_Classes_Download_Version: 1.4
# ELM_Classes_Download_Date: 2016-05-09 10:55:27.098737
# Origin: elm.eu.org
# Type: tsv
# Num_Classes: 247
ELMdefinitions = pd.read_csv(absolute_path + "/support/elm_classes.csv")


# PepX database motives
#
PEPXdefinitions = pd.read_csv(absolute_path + "/support/elms_pepx.csv:)


def motif_search(seq, ELMdefinitions=ELMdefinitions):
    """.

    Return all the ELM in a sequences.

    Parameters
    ----------
    seq: str
        sequence to query

    ELMdefinitions: pandas.DataFrame (default: totalELMdefinitions)
        Defintions and information about ELM in pandas format. Alternative
        PepXELMdefinitions can be used, ELM of domains in pepX

    Returns
    -------
        Returns a Datafram with all the motif in the sequnece
        columns {Accession,ELMIdentifier,Description,Regex,Probability,start,end}

    """
    elmhits = list()
    for i, elm in ELMdefinitions.iterrows():
        # regex = Def[i][1]
        m = re.search(elm['Regex'], seq)
        if str(m) != "None":
            elm['start'] = m.start()
            elm['end'] = m.end()
            elmhits.append(elm)

    return pd.DataFrame(elmhits)

# Test
# motif('A1L3X0', 'MNSVGEACTDMKREYDQCFNRWFAEKFLKGDSSGDPCTDLFKRYQQCVQKAIKEKEIPIEGLEFMGHGKEKPENSS')
# a = motifsearch('MNSVGEACTDMKREYDQCFNRWFAEKFLKGDSSGDPCTDLFKRYQQCVQKAIKEKEIPIEGLEFMGHGKEKPENSS\
    # PPPPPPPP')
# print a
# print a.shape

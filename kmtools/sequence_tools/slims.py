"""
Tools to find Short Linear Motifs on Protein sequences - SLiM -

"""
import os
import re
import logging

import pandas as pd

logger = logging.getLogger(__name__)
<<<<<<< HEAD
# To improve. Motif definitios repositories
absolute_path = os.path.dirname(__file__)


=======

# TODO: improve this
#Load Motif definitios repositories
absolute_path = os.path.dirname(__file__)
# File details:
>>>>>>> 38d96f4b411cd430164a4af927f4b4497f60f5ad
# ELM_Classes_Download_Version: 1.4
# ELM_Classes_Download_Date: 2016-05-09 10:55:27.098737
# Origin: elm.eu.org
# Type: tsv
# Num_Classes: 247
<<<<<<< HEAD
ELMdefinitions = pd.read_csv(absolute_path + "/support/elm_classes.csv")


# PepX database motives
#
PEPXdefinitions = pd.read_csv(absolute_path + "/support/elms_pepx.csv:)


def search_motif(seq, ELMdefinitions=ELMdefinitions):
=======
ELM_SLiM_defs = pd.read_csv(absolute_path + "/support/elm_classes.tsv", sep='\t',
                                  encoding='utf-8')
# Load
# PepX database motives
PepXELM_SLiM_defs = pd.read_csv(absolute_path + "/support/Pepx_completeELM_Redex.tsv", sep='\t')


def search_motif(seq, ELMdefinitions=ELM_SLiM_defs):
>>>>>>> 38d96f4b411cd430164a4af927f4b4497f60f5ad
    """Return all the Motif in a sequences


    Parameters
    ----------
    seq: str
        sequence to query

    ELMdefinitions: str, list, dict, pandas.Dataframe['Regex']
        Motif definitions, by default use  /support/elm_classes.tsv loaded to a pandas.DataFrame.
        Alternative you can use PepXELM_SLiM_defs - SLiM of domains in PepX.
        Custom regular expression can be inserted as str r'P.L', list [r'P.L',r'P.L.Y',r'L...P.L']
        or  dict {'PxL SHORT':r'P.L','Experimental PxLxY':r'P.L.Y'}


    Return
    ------
        Returns a Datafram with all the motif in the sequnece
        columns {Accession,Identifier,Description,Regex,Probability,start,end}

    """
    # Check input Parameter type and transform to pandas DataFrame
    if isinstance(ELMdefinitions, str):
        ELMdefinitions = list(ELMdefinitions)

    if isinstance(ELMdefinitions, list):
        ELMdefinitions = _convert_list(ELMdefinitions)

    if isinstance(ELMdefinitions, dict):
        ELMdefinitions = _convert_dict(ELMdefinitions)

    elmhits = list()
    for i, elm in ELMdefinitions.iterrows():
        # regex = Def[i][1]
        m = re.search(elm['Regex'], seq)
        if str(m) != "None":
            elm['start'] = m.start()
            elm['end'] = m.end()
            elmhits.append(elm)

    return pd.DataFrame(elmhits)


def _convert_list(slim_list):
    """Convert a list to a DataFrame.

    Parameters
    ----------

    Returns
    -------

    """
    # [r'P.L',r'P.L.Y',r'L...P.L']
    paried_list = [[i, x] for i, x in enumerate(slim_list)]
    return pd.DataFrame(paried_list, columns=['Identifier', 'Regex'])


def _convert_dict(slim_dict):
    """Convert a dictionary to a pandas DataFrame.

    Parameters
    ----------

    Returns
    -------

    """
    # {'PxL SHORT':r'P.L','Experimental PxLxY':r'P.L.Y'}
    return pd.DataFrame.from_dict(slim_dict,
                                  orient='index').reset_index().rename(columns={0: 'Regex',
                                                                                'index': 'Identifier'})


def _check_elm_dataframe(df):
    """Check DataFrame Format.

    Parameters
    ----------

    Returns
    -------

    """

    if 'Regex' in df.columns and df.shape[0] > 0:
        return True
    else:
        return False
<<<<<<< HEAD
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
=======
>>>>>>> 38d96f4b411cd430164a4af927f4b4497f60f5ad

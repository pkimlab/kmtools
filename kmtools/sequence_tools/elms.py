import os
import re
import pandas as pd


# To improve. Motif definitios repositories
absolute_path = os.path.dirname(__file__)


#ELM_Classes_Download_Version: 1.4
#ELM_Classes_Download_Date: 2016-05-09 10:55:27.098737
#Origin: elm.eu.org
#Type: tsv
#Num_Classes: 247
totalELMdefinitions = pd.read_csv(absolute_path+"/support/elm_classes.tsv", sep='\t',encoding='utf-8')


# PepX database motives
#
PepXELMdefinitions = pd.read_csv(absolute_path+"/support/Pepx_completeELM_Redex.tsv", sep='\t')




def motif_search(seq,ELMdefinitions=totalELMdefinitions):

    '''Return all the ELM in a sequences
    Parameters
    ----------
    seq: str
        sequence to query

    ELMdefinitions: pandas.DataFrame (default: totalELMdefinitions)
        Defintions and information about ELM in pandas format. Alternative
        PepXELMdefinitions can be used, ELM of domains in pepX

    Return
    ------
        Returns a Datafram with all the motif in the sequnece
        columns {Accession,ELMIdentifier,Description,Regex,Probability,start,end}

    '''

    elmhits = list()
    for i,elm in ELMdefinitions.iterrows():
        # regex = Def[i][1]
        m = re.search(elm['Regex'], seq)
        if str(m) != "None":
            elm['start'] = m.start()
            elm['end'] = m.end()
            elmhits.append(elm)

    return pd.DataFrame(elmhits)

# Test
# motif('A1L3X0', 'MNSVGEACTDMKREYDQCFNRWFAEKFLKGDSSGDPCTDLFKRYQQCVQKAIKEKEIPIEGLEFMGHGKEKPENSS')
# a = motifsearch('MNSVGEACTDMKREYDQCFNRWFAEKFLKGDSSGDPCTDLFKRYQQCVQKAIKEKEIPIEGLEFMGHGKEKPENSSPPPPPPPP')
# print a
# print a.shape




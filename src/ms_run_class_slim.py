
from pathlib import Path
import subprocess
import pandas as pd
import re
import numpy as np

CONTAMINANTS = ['ALBU_BOVIN',
 'AMYS_HUMAN',
 'CAS1_BOVIN',
 'CAS2_BOVIN',
 'CASB_BOVIN',
 'CASK_BOVIN',
 'CTRA_BOVIN',
 'CTRB_BOVIN',
 'K1C10_HUMAN',
 'K1C15_SHEEP',
 'K1C9_HUMAN',
 'K1H1_HUMAN',
 'K1H2_HUMAN',
 'K1H4_HUMAN',
 'K1H5_HUMAN',
 'K1H6_HUMAN',
 'K1H7_HUMAN',
 'K1H8_HUMAN',
 'K1HA_HUMAN',
 'K1HB_HUMAN',
 'K1M1_SHEEP',
 'K1M2_SHEEP',
 'K22E_HUMAN',
 'K2C1_HUMAN',
 'K2M1_SHEEP',
 'K2M2_SHEEP',
 'K2M3_SHEEP',
 'KRA33_SHEEP',
 'KRA34_SHEEP',
 'KRA3A_SHEEP',
 'KRA3_SHEEP',
 'KRA61_SHEEP',
 'KRB2A_SHEEP',
 'KRB2B_SHEEP',
 'KRB2C_SHEEP',
 'KRB2D_SHEEP',
 'KRHB1_HUMAN',
 'KRHB2_HUMAN',
 'KRHB3_HUMAN',
 'KRHB4_HUMAN',
 'KRHB5_HUMAN',
 'KRHB6_HUMAN',
 'KRUC_SHEEP',
 'LALBA_BOVIN',
 'LYSC_LYSEN',
 'PEPA_BOVIN',
 'PEPA_PIG',
 'PEPB_PIG',
 'PEPC_PIG',
 'SSPA_STAAU',
 'TRY1_BOVIN',
 'TRY2_BOVIN',
 'TRYP_PIG',
 'ALBU_HUMAN',
 'ANT3_HUMAN',
 'ANXA5_HUMAN',
 'B2MG_HUMAN',
 'BID_HUMAN',
 'CAH1_HUMAN',
 'CAH2_HUMAN',
 'CATA_HUMAN',
 'CATD_HUMAN',
 'CATG_HUMAN',
 'CO5_HUMAN',
 'CRP_HUMAN',
 'CYB5_HUMAN',
 'CYC_HUMAN',
 'EGF_HUMAN',
 'FABPH_HUMAN',
 'GELS_HUMAN',
 'GSTA1_HUMAN',
 'GSTP1_HUMAN',
 'HBA_HUMAN',
 'HBB_HUMAN',
 'IGF2_HUMAN',
 'IL8_HUMAN',
 'KCRM_HUMAN',
 'LALBA_HUMAN',
 'LEP_HUMAN',
 'LYSC_HUMAN',
 'MYG_HUMAN',
 'NEDD8_HUMAN',
 'NQO1_HUMAN',
 'NQO2_HUMAN',
 'PDGFB_HUMAN',
 'PPIA_HUMAN',
 'PRDX1_HUMAN',
 'RASH_HUMAN',
 'RETBP_HUMAN',
 'SODC_HUMAN',
 'SUMO1_HUMAN',
 'SYH_HUMAN',
 'TAU_HUMAN',
 'THIO_HUMAN',
 'TNFA_HUMAN',
 'TRFE_HUMAN',
 'TRFL_HUMAN',
 'UB2E1_HUMAN',
 'UBE2C_HUMAN',
 'UBE2I_HUMAN',
 'RS27A_HUMAN',
 'GAG_SCVLA',
 'CYC_HORSE',
 'CAH2_BOVIN',
 'ADH1_YEAST',
 'ALDOA_RABIT',
 'LYSC_CHICK',
 'MYG_HORSE',
 'OVAL_CHICK',
 'BGAL_ECOLI',
 'DHE3_BOVIN',
 'GFP_AEQVI',
 'SRPP_HEVBR',
 'REF_HEVBR',
 'PLMP_GRIFR',
 'KKA1_ECOLX']


def setClass(accession, contaminants_list=CONTAMINANTS):
    '''
    assign a class to the peptide
    NEOTRANSCRIPT is the tag that is added to neotranscript entries by the add_PE.py script
    lcl is a tag added to the result of translation with ORFfinder
    sp is the tag of proteins from swissprot (uniprot)
    source of the contaminants list: https://www.thegpm.org/crap/
    '''
    if ((('NEOTRANSCRIPT' in accession) or ('lcl' in accession)) and 'sp' not in accession):
        return("Non canonical")
    elif ((('NEOTRANSCRIPT' in accession) or ('lcl' in accession)) and 'sp' in accession):
        return("Ambiguous")
    elif any(contaminant in accession for contaminant in contaminants_list):
        return('Contaminant')
    else:
        return('Canonical')


class MSrun:
    def __init__(self, metadata_row):
        self.sample = metadata_row["sample"]
        self.DB = metadata_row["DB"]
        self.search_engine = metadata_row["search_engine"]
        self.fdr = metadata_row["fdr"]
        if "hla_alleles" in metadata_row:
            self.hla = metadata_row["hla_alleles"].split("|")

        search_engine = metadata_row["search_engine"]
        print(metadata_row["path"])
        if "Fragpipe" in search_engine:
            self.psm = self.from_fragpipe(metadata_row["path"])
        else:
            raise ValueError(f"Provided search engine ({search_engine}) is not supported.\n \
            Only Fragpipe search engine is supported.")
        self.unique = self.psm.drop_duplicates(subset=['seq_clear'], inplace=False)

    @staticmethod
    def from_fragpipe(path):
        def get_fraction(row):
            return row["Spectrum"].split(".")[0]
	
        def get_scan(row):
            return int(row["Spectrum"].split(".")[2])

        def extract_accessions(row):
            if not (row["Is Unique"] | isinstance(row["Mapped Proteins"], float)):
                accessions = ",".join((row["Protein"],row["Mapped Proteins"]))
                accessions = accessions.replace(",",";")
                accessions = accessions.replace(" ","")
            else:
                accessions = row["Protein"]
            return str(accessions)

        psm_file = pd.read_csv(path, sep="\t")
        #print(len(psm_file.index))
        dict_info = {"type": ["PEPTIDE"]*len(psm_file),
            "rt": psm_file["Retention"], # convert into seconds
            #"expectation": psm_file["Expectation"], # see http://pappso.inrae.fr/bioinfo/i2masschroq/documentation/html/chap_fundamentals-bottom-up-proteomics.html E-value
            #"spectral_entropy": psm_file["SpectralSim"], # another quality metric for spectra see https://doi.org/10.1038/s41592-021-01331-z 
            "score": psm_file["Hyperscore"], # using PEP score as q-value
            "rank": [None]*len(psm_file),
            "sequence": psm_file["Modified Peptide"],
            "seq_clear": psm_file["Peptide"],
            "aa_before": psm_file["Prev AA"],
            "aa_after": psm_file["Next AA"],
            "score_type": ["PSM"]*len(psm_file),
            "search_identifier": ["Fragpipe"]*len(psm_file),
            "accessions": psm_file.apply(extract_accessions, axis=1),
            "start": psm_file["Protein Start"],
            "end": psm_file["Protein End"],
            "fraction": psm_file.apply(get_fraction, axis=1),
            "scan": psm_file.apply(get_scan, axis=1),
            #"intensities": psm_file["Intensities"]
            #    "sample": meta["sample"],
            "fdr": [None]*len(psm_file)
            }
        psm = pd.DataFrame.from_dict(dict_info)
        return psm

    def isCanonical(self):
        self.psm["class"] = self.psm.accessions.apply(lambda x: setClass(x))
        self.unique["class"] = self.unique.accessions.apply(lambda x: setClass(x))
import pandas as pd
from Bio import SeqIO
import ahocorasick
import re
import itertools
import argparse
import time

def clear_seq(s):
    s_2 = re.sub(r"[^A-Z]", "", s)
    return s_2

def mutate_L_to_I_and_I_to_L(seq, pos_to_change):
    '''
    changes all I and L in the string to the one another at pos_to_change
    '''
    for pos in pos_to_change:
        if seq[pos] == 'I':
            seq = seq[:pos] + 'L' + seq[pos+1:]
        elif seq[pos] == 'L':
            seq = seq[:pos] + 'I' + seq[pos+1:]
        else:
            print(f'Given string does not equal I or L at position {pos}, {seq[pos]} was found instead!')
    return seq

def all_combinations(any_list):
    return itertools.chain.from_iterable(
        itertools.combinations(any_list, i + 1)
        for i in range(len(any_list)))

def get_all_IL_LI_substitutions(peptide):
    mutated_peptides = []
    IL_pos = [m.start() for m in re.finditer(r'[IL]', peptide)]
    IL_pos_all_comb = [list(l) for l in all_combinations(IL_pos)]
    for comb in IL_pos_all_comb:
        peptide_mutated = mutate_L_to_I_and_I_to_L(peptide, comb)
        mutated_peptides.append(peptide_mutated)
    return mutated_peptides

def search_in_database_with_I_L_changed(peptides, database):
    """
    Match each peptide with database entries in uniprot using PyAhoCorasick.
    
    Args:
        peptides (list): List of peptides to search for.
        uniprot (or other protein fasta) (DataFrame): DataFrame with 'seq' (sequences) and 'name' columns.
        
    Returns:
        dict: A dictionary where keys are peptides and values are lists of matching database entries.
    Also searches all possible peptides derived from the given list by changing I to L and L to I
    """
    automaton = ahocorasick.Automaton()
    mutated_peptides = []
    for peptide in peptides:
        automaton.add_word(peptide, peptide)
        mutated_peptides.extend(get_all_IL_LI_substitutions(peptide))
    for peptide in mutated_peptides:
        automaton.add_word(peptide, peptide)
    peptides.extend(mutated_peptides)
    automaton.make_automaton()
    results = {peptide: [] for peptide in peptides}

    for index, row in database.iterrows():
        sequence = str(row['seq'].seq)  
        db_name = row['name']         


        for end_idx, peptide in automaton.iter(sequence):
            if peptide not in results:
                results[peptide] = []
            results[peptide].append(db_name)

    return results

def main():
    parser = argparse.ArgumentParser(description='Map Casanovo ouput peptides to fasta database(s)')
    parser.add_argument('--input_mztab', required=True, help='Path to Casanovo output .mztab')
    parser.add_argument('--output_csv', required=True, help='Path to output .csv')
    parser.add_argument('--canonical_proteome', required=True, help='Path to canonical proteome fasta')
    parser.add_argument('--extra_fasta', required=False, help='Path to extra fasta to include in search')
    args = parser.parse_args()
    #0 open the canonical proteome as dataframe
    uniprot = SeqIO.to_dict(SeqIO.parse(args.canonical_proteome, "fasta"))
    d = {'name': uniprot.keys(), 'seq': uniprot.values()}
    df_uni = pd.DataFrame.from_dict(d)
    #1 open the input mztab as dataframe
    with open(args.input_mztab,"r") as fi:
        k=0
        for ln in fi:
            k+=1
            if ln.startswith("PSH"):
                starts=k-1
    df = pd.read_csv(args.input_mztab, sep='\t', header=starts)
    #2 map peptides to the canonical proteome
    start_time = time.time()
    df['seq_clear'] = df['sequence'].apply(clear_seq)
    search_result = search_in_database_with_I_L_changed(df['seq_clear'].unique().tolist(), df_uni)
    df['IL_substitution'] = df['seq_clear'].apply(get_all_IL_LI_substitutions)
    df = df.explode('IL_substitution')
    df['canonical_proteome_match'] = df['seq_clear'].apply(lambda x: search_result[x])
    df['canonical_proteome_match_with_IL_change'] = df['IL_substitution'].apply(lambda x: search_result[x] if not(pd.isna(x)) else [])
    df['canonical_proteome_match_bool'] = df['canonical_proteome_match'].apply(lambda x: True if len(x)!=0 else False)
    df['canonical_proteome_bool_with_IL_change'] = df['canonical_proteome_match_with_IL_change'].apply(lambda x: True if len(x)!=0 else False)
    #3 if extra fasta is given, map to that too
    if args.extra_fasta:
        extra_db = SeqIO.to_dict(SeqIO.parse(args.extra_fasta, "fasta"))
        d_extra = {'name': extra_db.keys(), 'seq': extra_db.values()}
        df_extra = pd.DataFrame.from_dict(d_extra)
        search_result_extra = search_in_database_with_I_L_changed(df['seq_clear'].unique().tolist(), df_extra)
        df['extra_fasta_match'] = df['seq_clear'].apply(lambda x: search_result_extra[x])
        df['extra_fasta_match_with_IL_change'] = df['IL_substitution'].apply(lambda x: search_result_extra[x] if not(pd.isna(x)) else [])
        df['extra_fasta_match_bool'] = df['extra_fasta_match'].apply(lambda x: True if len(x)!=0 else False)
        df['extra_fasta_bool_with_IL_change'] = df['extra_fasta_match_with_IL_change'].apply(lambda x: True if len(x)!=0 else False)
    df.to_csv(args.output_csv, index=False)
    print("Process finished --- %s seconds ---" % (time.time() - start_time))

if __name__=="__main__":
    main()

import sys
import os
import numpy as np
import pandas as pd
import re

from progress import Progress


def read_swaths2stats_peptide_table( event, 
                                    log_fh, 
                                    progressData,
                                    in_filename):

    progress = Progress(progressData, "percentage-indicator")

    step = 0
    n_steps = 2
    df = pd.read_csv(in_filename, sep='\t')
    step += 1
    progress.update_n_of_m(step, n_steps)
    df["Peptide"] = df["ProteinName_FullPeptideName"].apply(lambda x: x.split("_")[-1].strip())
    df["Proteins"] = df["ProteinName_FullPeptideName"].apply(lambda x: ";".join(set("_".join(x.split("_")[0:-1]).split("/")[1:])) )
    step += 1
    progress.update_n_of_m(step, n_steps)
    progress.ready()

    return df

def annotate_df(
        event, 
        log_fh, 
        progressData,
        df,
        dicts,
        out_matrix_filename,
        out_annot_filename,
        ambiguous_threshold,
        merge_unimods,
        contaminants):

    progress = Progress(progressData, "percentage-indicator")

    sorted_samplenames = df.columns.tolist()
    sorted_samplenames.remove('Proteins')
    sorted_samplenames.remove('Peptide')
    sorted_samplenames.remove("ProteinName_FullPeptideName")

    if merge_unimods:
        print ("Unimplemented")

        # rules_dict = {"Proteins" : lambda x: ";".join(x), 
        #               "FullPeptideName" : lambda x: ";".join(x),}
        # for filename in sorted_samplenames:
        #                 rules_dict[filename] = "sum"

        # df = df.groupby("Sequence", as_index=False).agg(rules_dict)
        # df['Proteins'] = df['Proteins'].apply(lambda x: ";".join(set(x.split(';'))))
        # df['FullPeptideName'] = df['FullPeptideName'].apply(lambda x: ";".join(set(x.split(';'))))

    nrows, _ = df.shape

    step = 0
    n_steps = nrows

    annot_colnames = list(dicts)

    for annot_colname in annot_colnames:
        df[annot_colname]=["" for x in range(nrows)]


    print ("Iterating rows")
    for index, row in df.iterrows():

        sets = {}
        
        for colname in annot_colnames:
            if colname not in sets:
                sets[colname] = set()

        for ID in row["Proteins"].split(";"):
            for k in dicts:
                if ID in dicts[k]:
                    sets[k].add(dicts[k][ID])

        L = {}
        for k in annot_colnames:
            L[k] = list(sets[k])

            if len(L[k]) > 1 and "unknown" in L[k]:
                L[k].remove('unknown')

            if len(L[k]) == 1:
                df.at[df.index[index], k] = L[k][0]
            elif len(L[k]) > 1:
                if not ambiguous_threshold or (ambiguous_threshold and len(L[k]) < int(ambiguous_threshold)):
                    df.at[df.index[index], k] = ";".join(L[k])
                else:
                    df.at[df.index[index], k] = "ambiguous"
            else:
                df.at[df.index[index], k] = "unknown"

        step += 1
        progress.update_n_of_m(step, n_steps)
        if event and event.kill_flag:
            return

    if contaminants:
        df = df[~df.Proteins.str.contains("|".join(contaminants))]


    if out_annot_filename:
        
        print ("making annot file")

        annot_df_cols = ["ProteinName_FullPeptideName",'Peptide','Proteins']
        annot_df_cols.extend(annot_colnames)
        annot_df = df.filter(items=annot_df_cols)
        #annot_df.rename(columns={use_col:"FullPeptideName"}, inplace=True)
        annot_df.to_csv(out_annot_filename, sep="\t", index=False)

    if out_matrix_filename:
        print ("making matrix file")

        df_cols = ["ProteinName_FullPeptideName"] 
        df_cols.extend(sorted_samplenames)
        df = df.filter(df_cols)
        #df.rename(columns={use_col:"FullPeptideName"}, inplace=True)
        df.to_csv(out_matrix_filename, sep="\t", index=False)

    progress.update_n_of_m(n_steps, n_steps)
    progress.ready()

    return 

#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import os
import sys
import argparse

def readPDB(PDB_ID):
    def pdbParseClean(line):

        if line.startswith("ATOM"):
            field_record = "ATOM"
        if line.startswith("HETATM"):
            field_record = "HETATM"

        field_chain   = line[21]
        field_serial  = line[6:12]
        field_notsure = line[12]
        field_atom    = line[13:16]
        field_altLoc  = line[16]
        field_resname = line[17:20]

        field_inser   = line[26]
        field_resseq  = line[23:26]
        field_x       = line[30:38]
        field_y       = line[38:46]
        field_z       = line[46:54]
        field_occ     = line[54:60]
        field_temp    = line[60:67]
        field_element = line[76:80]

        return [field_record,field_serial,field_atom,field_altLoc,field_resname,field_chain,field_resseq,field_inser,field_x,field_y,field_z,field_occ,field_temp,field_element,field_notsure]

    pdb_array = []

    pdb_columns = ["record","serial","atom","altLoc","resname","chain","resseq","insertion","x","y","z","occ","temp","element","notsure"]
    
    with open(PDB_ID,"r") as infFile:
        for line in infFile:
            line = line.replace("\n","")

            if line.startswith("ATOM") or line.startswith("HETATM"):

                line_array = pdbParseClean(line)
                line_array = [i.replace(" ","") for i in line_array]

                if len(line) > 80:
                        field_charge = line[78:80]
                        line_array.append(field_charge)
                        pdb_columns.append("charge")

                pdb_array.append(line_array)


    pdb_df = pd.DataFrame.from_records(pdb_array, columns=pdb_columns)

    return pdb_df
    
def writePDB(df,name):
    def writelines(row,cleanChain):    
        write_string = row[0].ljust(6)+row[1].rjust(5)+" "

        if row[14] != "":
            write_string = write_string+row[14]
        else:
            write_string = write_string+" "

        write_string = write_string+row[2].ljust(3)+row[3].ljust(1)+row[4].ljust(4)+row[5].ljust(2)+row[6].rjust(3)+row[7].ljust(1)+"   "+row[8].rjust(8)+row[9].rjust(8)+row[10].rjust(8)+"  "+row[11].ljust(3)+" "+row[12].rjust(5)

        if len(row[12]) > 5:
            write_string = write_string+"          "
        else:
            write_string = write_string+"           "

        write_string=write_string+row[13].ljust(3)+"\n"

        cleanChain.write(write_string)

    pdbName      = name+".pdb" 
    outputPDB    = open(pdbName,"w")
    df["serial"] = df["serial"].astype(str)

    for row in df.itertuples(index=False):
        writelines(row,outputPDB)

    outputPDB.close()

def individual_excipients(df):
    unique_exp = df.resseq.unique()

    unique_df_array = []
    unique_id_array = []

    for resseqid in unique_exp:
        unique_df = df[df["resseq"] == resseqid].copy(deep=True)
        unique_df_array.append(unique_df)
        unique_id_array.append(resseqid)

    return unique_df_array, unique_id_array


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--PDB", type=str, help="Input file name")
    parser.add_argument("--ligands", type=str, help="ligands to be separated")
    args = parser.parse_args()

    inputfile             = args.PDB
    excipient_list        = args.ligands.split(",")
    unique_excipient_list = []

    PDB_basename   = os.path.basename(inputfile)
    PDB_ID         = os.path.splitext(PDB_basename)[0]

    if os.path.isfile(inputfile):
        try:
            whole_pdb_df = readPDB(inputfile)
        except:
            print(f"Couldn't read {PDB_ID}")
        else:
            protein_df     =  whole_pdb_df[whole_pdb_df["record"] == "ATOM"].copy(deep=True)

            if len(protein_df) < 100000:
                writePDB(protein_df,PDB_ID+"_protein")
                HET_df     =  whole_pdb_df[whole_pdb_df["record"] == "HETATM"].copy(deep=True)

                print(f"{PDB_ID}")
                print(f"Number of protein atoms: {len(protein_df)}")

                for excipient in excipient_list:
                    exp_df = HET_df[HET_df["resname"] == excipient].copy(deep=True)

                    unique_df_array, unique_id_array = individual_excipients(exp_df)

                    if len(unique_id_array) == 0:
                        print(f"No excipient {excipient} found")
                    else:
                        print(f"Excipient {excipient} found")

                        for idx, df in enumerate(unique_df_array):
                            combined_df = pd.concat([protein_df,df])
                            writePDB(combined_df,PDB_ID+"_protein_excipient_"+excipient+"_"+unique_id_array[idx])
                            unique_excipient_list.append(excipient+"_"+unique_id_array[idx])

                print(f"Excipients: {' '.join(str(x) for x in unique_excipient_list)}\n")

            else:
                print("PDB too large, not writing")
    else:
        print(f"File {inputfile} not found")

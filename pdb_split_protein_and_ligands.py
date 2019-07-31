#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import os
import sys
import argparse

def readPDB(PDB_ID, parse_type):
    
    def pdbParseClean(line,atom_type):
        field_chain = line[21].replace(" ","")
        field_record = atom_type
        field_serial = line[6:12].replace(" ","")
        field_notsure = line[12].replace(" ","")
        field_atom = line[13:16].replace(" ","")
        field_altLoc = line[16].replace(" ","")
        field_resname = line[17:20].replace(" ","")

        field_inser = line[26].replace(" ","")
        field_resseq = line[23:26].replace(" ","")
        field_x = float(line[30:38].replace(" ",""))
        field_y = float(line[38:46].replace(" ",""))
        field_z = float(line[46:54].replace(" ",""))
        field_occ = line[54:60].replace(" ","")
        field_temp = line[60:67].replace(" ","")
        field_element = line[76:80].replace(" ","")

        return field_record,field_serial,field_atom,field_altLoc,field_resname,field_chain,field_resseq,field_inser,field_x,field_y,field_z,field_occ,field_temp,field_element,field_notsure

    
    def parsePDB(PDB_ID,atom_type):
        
        pdb_name = PDB_ID+".pdb"
        pdb = []
        
        with open(pdb_name,"r") as infFile:
            for line in infFile:
                line = line.replace("\n","")
                if line.startswith(atom_type):
                    field_record,field_serial,field_atom,field_altLoc,field_resname,field_chain,field_resseq,field_inser,field_x,field_y,field_z,field_occ,field_temp,field_element,field_notsure = pdbParseClean(line, atom_type)

                    field_array = [field_record,field_serial,field_atom,field_altLoc,field_resname,field_chain,field_resseq,field_inser,field_x,field_y,field_z,field_occ,field_temp,field_element,field_notsure]
                    pdb_columns = ["record","serial","atom","altLoc","resname","chain","resseq","insertion","x","y","z","occ","temp","element","notsure"]

                    if len(line) > 80:
                            field_charge = line[78:80]
                            field_array.append(field_charge)
                            pdb_columns.append("charge")

                    pdb.append(field_array)


        pdb_df = pd.DataFrame.from_records(np.array(pdb), columns=pdb_columns)

        return pdb_df
    
    if parse_type == "protein":
        return parsePDB(PDB_ID,"ATOM")
    if parse_type == "excipient":
        return parsePDB(PDB_ID,"HETATM")


def writePDB(df,name):
    def writelines(row,cleanChain):    
        #cleanChain.write(row)
        crow = row.str.strip()
        #ATOM record
        cleanChain.write(crow[0].ljust(6))
        #Serial
        cleanChain.write(crow[1].rjust(5))
        cleanChain.write(" ")
        if crow[14] != "":
            cleanChain.write(crow[14])
        else:
            cleanChain.write(" ")
        #Atom name
        cleanChain.write(crow[2].ljust(3))
        #Altloc
        cleanChain.write(crow[3].ljust(1))
        #resname
        cleanChain.write(crow[4].ljust(4))
        #chain
        cleanChain.write(crow[5].ljust(2))
        #number
        cleanChain.write(crow[6].rjust(3))
        #if any(c.isalpha() for c in crow[6]):
        cleanChain.write("")

        cleanChain.write(crow[7].ljust(1))
        cleanChain.write("   ")

        #xyz
        cleanChain.write(crow[8].rjust(8))
        cleanChain.write("")
        cleanChain.write(crow[9].rjust(8))
        cleanChain.write("")
        cleanChain.write(crow[10].rjust(8))
        cleanChain.write("  ")

        #occ
        cleanChain.write(crow[11].ljust(3))
        cleanChain.write(" ")
        #temp
        cleanChain.write(crow[12].rjust(5))
        if len(crow[12]) > 5:
            cleanChain.write("          ")
        else:
            cleanChain.write("           ")
        #element
        cleanChain.write(crow[13].ljust(3))
        cleanChain.write("\n")

    pdbName      = name+".pdb" 
    outputPDB    = open(pdbName,"w")
    df["serial"] = df["serial"].astype(str)

    for i, row in df.iterrows():
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
            HET_df = readPDB(PDB_ID,parse_type="excipient")
        except:
            print(f"{PDB_ID} no ligands in pdb")
        else:
            protein_df = readPDB(PDB_ID,parse_type="protein")
            writePDB(protein_df,PDB_ID+"_protein")

            for excipient in excipient_list:
                exp_df = HET_df[HET_df["resname"] == excipient].copy(deep=True)

                if len(exp_df) > 0:
                    unique_df_array, unique_id_array = individual_excipients(exp_df)

                    for idx, df in enumerate(unique_df_array):
                        writePDB(df,PDB_ID+"_excipient_"+excipient+"_"+unique_id_array[idx])
                        combined_df = pd.concat([protein_df,df])
                        writePDB(combined_df,PDB_ID+"_protein_excipient_"+excipient+"_"+unique_id_array[idx])
                        unique_excipient_list.append(excipient+"_"+unique_id_array[idx])

                        print(f"Making pdb files for excipient {excipient}{unique_id_array[idx]}")
                        writePDB(exp_df,PDB_ID+"_excipient_"+excipient)
                else:
                    print(f"{excipient} not present in structure")
    else:
        print(f"File {inputfile} not found")

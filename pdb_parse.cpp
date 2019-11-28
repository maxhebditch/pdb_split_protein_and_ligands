#include <string>
#include <iostream>
#include <iomanip>
#include <bits/stdc++.h>
#include <stdio.h>
#include <algorithm>

using namespace std;

string strip_space(string str){
        str.erase(remove(str.begin(), str.end(), ' '), str.end());
        return str;
}

int main(int argc, char* argv[]) {

        /* check write number of args passed */
        if (argc < 1) {
                std::cerr << "Usage: " << argv[0] << " protein.sa" << std::endl;
                return 1;
        }

        /* stream the input file */
        string pdb_file = argv[1];
        ifstream pdb_stream (pdb_file);
        string pdb_line;

        /* initialise variables to hold ASA details */
        string record_name;
        int serial;
        string notsure;
        string AA_atom;
        string AA_altloc;
        string AA_resname;
        string AA_chain;
        int AA_resseq;
        string AA_insertion;
        float AA_x;
        float AA_y;
        float AA_z;
        float AA_occ;
        float AA_temp;
        string AA_element;
        string AA_charge;

        /* go through file line by line */
        while(getline(pdb_stream, pdb_line )) {
                /* if line starts with > */
                record_name = strip_space(pdb_line.substr(0,6));
                if (record_name == "ATOM"){

                        serial       = stoi(pdb_line.substr(6,5));
                        notsure      = pdb_line.substr(12,1);
                        AA_atom      = strip_space(pdb_line.substr(13,4));
                        AA_altloc    = pdb_line.substr(16,1);
                        AA_resname   = strip_space(pdb_line.substr(17,3));
                        AA_chain     = pdb_line.substr(21,1);
                        AA_resseq    = stoi(pdb_line.substr(22,4));
                        AA_insertion = pdb_line.substr(26,1);
                        AA_x         = stof(pdb_line.substr(30,8));
                        AA_y         = stof(pdb_line.substr(38,8));
                        AA_z         = stof(pdb_line.substr(46,8));
                        AA_occ       = stof(pdb_line.substr(54,6));
                        AA_temp      = stof(pdb_line.substr(60,6));
                        AA_element   = strip_space(pdb_line.substr(76,2));
                        AA_charge    = strip_space(pdb_line.substr(78,2));

                        /* col  1-6  RECORD NAME */
                        cout << left  << setw(6) << record_name;
                        /* col  7-11 SERIAL */
                        cout << right << setw(5) << serial; 
                        /* col    12 EMPTY */
                        cout << " ";
                        cout << notsure;
                        /* col 13-16 Atom name */
                        cout << left << setw(3) << AA_atom;
                        /* col    17 Altloc */
                        cout << AA_altloc;
                        /* col 18-20 resname */
                        cout << right << setw(3) << AA_resname;
                        /* col    21 EMPTY */
                        cout << " ";
                        /* col    22 CHAIN ID */
                        cout << AA_chain;
                        /* col 23-26 resseq */
                        cout << right << setw(4) << AA_resseq;
                        /* col    27 CHAIN ID */
                        cout << AA_insertion;
                        /* col 28-30 empty */
                        cout << "   ";
                        /* col 31-38 X coords */
                        cout << fixed << setprecision(3);
                        cout << right << setw(8) << AA_x;
                        /* col 39-46 Y coords */
                        cout << right << setw(8) << AA_y;
                        /* col 47-54 Z coords */
                        cout << right << setw(8) << AA_z;
                        /* col 55-60 occupency */
                        cout << fixed << setprecision(2);
                        cout << right << setw(6) << AA_occ;
                        /* col 61-66 tempfactor */
                        cout << right << setw(6) << AA_temp;
                        /* col 67-76 empty */
                        cout << "          ";
                        /* col 77-78 element */
                        cout << right << setw(2) << AA_element;
                        /* col 79-80 charge */
                        cout << right << setw(2) << AA_charge;

                        cout << endl;
                }
        }



        return 0;

}

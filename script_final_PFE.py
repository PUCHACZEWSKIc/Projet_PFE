from Extract_fasta import Fasta_extract
from gap_insertion import Gap_Insertion
import argparse
import glob
import subprocess
import os 
import numpy as np


if __name__ == "__main__" : 
    arg_manager = argparse.ArgumentParser()
    arg_manager.add_argument('-infile', '-i', type = str)
    arg_manager.add_argument('-add_new_sequences', '-n', type = str)
    arg_manager.add_argument('-pack_size', '-ps', type = int)
    args = arg_manager.parse_args()
    if args.add_new_sequences == "y": 
        file = Fasta_extract(args.infile)
        file.sequences()
        file.add_sequences()
    else : 
        file = Fasta_extract(args.infile)
        file.sequences()
        file.sequences_packages(args.pack_size)
    
    if not os.path.exists("sequences_packages_align"):
        os.makedirs("sequences_packages_align")
    if not os.path.exists("sequences_packages_compo"):
        os.makedirs("sequences_packages_compo")
    if not os.path.exists("sequences_packages_cons"):
        os.makedirs("sequences_packages_cons")

    list_fichier = glob.glob('sequences_packages/pack_*.fasta')
    for subfile in range(len(list_fichier)) : 
        command_halign = "halign -o sequences_packages_align/aligned_pack_" + str(subfile) + ".fasta sequences_packages/pack_" + str(subfile) + ".fasta"
        subprocess.run(command_halign, shell = True, executable = "/bin/bash")
        fasta_int_infile = "sequences_packages_align/aligned_pack_" + str(subfile) + ".fasta"
        command_cons = "python3 compo_calc.py " + fasta_int_infile + " sequences_packages_compo/compo_" + str(subfile) + ".csv sequences_packages_cons/cons_" + str(subfile) + ".fasta"
        subprocess.run(command_cons, shell = True, executable = "/bin/bash")
    subprocess.run("rm -rd sequences_packages", shell = True, executable = "/bin/bash")
    
    list_cons = glob.glob("sequences_packages_cons/cons_*.fasta")
    list_cons_str = []
    for cons in list_cons : 
        txt = ""
        with open(cons) as read : 
            for line in read : 
                txt += line
        list_cons_str.append(txt)
    with open("all_cons.fasta", "a") as fileOut : 
        for cons in list_cons_str : 
            fileOut.write(cons)

    command_halign = "halign -o consensus_align.fasta all_cons.fasta"
    subprocess.run(command_halign, shell = True, executable = "/bin/bash")
    subprocess.run("rm -rd sequences_packages_align", shell = True, executable = "/bin/bash")
    subprocess.run("rm -rd sequences_packages_cons", shell = True, executable = "/bin/bash")
    subprocess.run("rm all_cons.fasta", shell = True, executable = "/bin/bash")


    list_pack_compos = glob.glob("sequences_packages_compo/compo_*.csv")
    for pack_compos in list_pack_compos : 
        data = np.genfromtxt(pack_compos, delimiter=',')
        compt = 0
        list_compos = []
        compos = []
        for i in data : 
            compos.append(i)
            if compt - 5 == 0 : 
                list_compos.append(compos)
                compt = 0
                compos = []
            else : 
                compt += 1
    
    subprocess.run("rm -rd sequences_packages_compo", shell = True, executable = "/bin/bash")

    command_cons = "python3 compo_calc.py consensus_align.fasta compo_all_cons.csv cons_all_cons.fasta"
    subprocess.run(command_cons, shell = True, executable = "/bin/bash")
    subprocess.run("rm consensus_align.fasta", shell = True, executable = "/bin/bash")

    extract_cons_fin = Fasta_extract("cons_all_cons.fasta")
    extract_cons_fin.sequences()

    subprocess.run("rm compo_all_cons.csv", shell = True, executable = "/bin/bash")
    subprocess.run("rm cons_all_cons.fasta", shell = True, executable = "/bin/bash")

    list_compos_final = []
    for i in list_compos :
        gap_inserter = Gap_Insertion()
        list_compos_final.append(gap_inserter.align_oldcons_cons(extract_cons_fin.sequence_list[0], i))
    compo_final = gap_inserter.fusion(list_compos_final)[0]

    with open("compos_cons_final.csv", 'w') as file_out: 
        np.savetxt(file_out, compo_final, delimiter=',', fmt='%d') 

    COMPOSITION = ["A", "C", "G", "T", "N", "-"]
    fin_cons = ''
    for i in compo_final.T : 
        fin_cons += [COMPOSITION[i.argmax()]][0]


    with open("cons_final.fasta", "a") as file_out : 
        file_out.write(fin_cons)

    





    


# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: chemtoxicity
#     language: python
#     name: chemtoxicity
# ---

import os
import sys
import glob
import pandas as pd

# +
#Get list of all pdb files
all_pdb_files = sorted(glob.glob("HumanProt_AF2_domains/*.pdb"))
all_pdb_files

for i in range(0,len(all_pdb_files)):
    if (i%1000==0):
        print("Completed "+str(i)+" proteins")
    
    pdb_file = all_pdb_files[i]
    #Convert one pdb file to pdbqt
    sample_name = pdb_file.split("/")[1]
    rev_sample_name_part1 = "-".join(sample_name.split("-")[0:2])
    rev_sample_name_part2 = "-".join(sample_name.split("_")[-2:])
    final_sample_name = rev_sample_name_part1+"-"+rev_sample_name_part2[:-4]
    print(final_sample_name)

    #Prepare the receptor
    command1 = "prepare_receptor -r "+pdb_file+" -A checkhydrogens -o "+final_sample_name+".pdbqt"
    os.system(command1)

    #Prepare the autosite
    command2 = "autosite -r "+final_sample_name+".pdbqt -o temp"
    os.system(command2)

    #Move rigidreceptor to alphafold id in pdbqt
    command3 = "mv temp/rigidReceptor.pdbqt "+final_sample_name+".pdbqt"
    os.system(command3)

    #Get the xyz corodinates for rigidReceptor
    fp = open("temp/rigidReceptor.gpf","r")
    all_lines = fp.readlines()
    xyz_coordinates = all_lines[8].strip("\r\n").split(" ")[1:]
    print(xyz_coordinates)
    fp.close()

    #Remove temporary files from AutoSite
    command4 = "rm -rf temp/* temp_AutoSiteSummary.log"
    os.system(command4)

    #If a folder not there for sample then create one
    if not os.path.exists(final_sample_name):
        os.makedirs(final_sample_name)

    #Move the pdbqt file to sample folder
    command5 = "mv "+final_sample_name+".pdbqt "+final_sample_name+"/"
    os.system(command5)

    #Create a drug folder exists in sample folder 
    drug_folder_name = final_sample_name+"/"+final_sample_name+"_drug_pdbqt"
    if not os.path.exists(drug_folder_name):
        os.makedirs(drug_folder_name)

    #Copy optimized ligands there
    command6 = "cp drug_pdbqt/* "+drug_folder_name
    os.system(command6)

    #Make a copy of sample configuration in the sample folder
    command7 = "cp sample_config.txt "+final_sample_name+"/"+final_sample_name+"_config.txt"
    os.system(command7)

    #Open the config file and make changes
    fp = open(final_sample_name+"/"+final_sample_name+"_config.txt","r")
    all_lines = fp.readlines()
    fp.close()

    pwd = os.getcwd()
    all_lines[0] = 'receptor = '+pwd+'/'+final_sample_name+"/"+final_sample_name+".pdbqt\n"
    all_lines[1] = 'ligand_directory = '+pwd+"/"+drug_folder_name+"\n"
    all_lines[2] = 'opencl_binary_path = '+os.path.dirname(pwd)+" \n"
    all_lines[3] = 'center_x = '+xyz_coordinates[0]+"\n"
    all_lines[4] = 'center_y = '+xyz_coordinates[1]+"\n"
    all_lines[5] = 'center_z = '+xyz_coordinates[2]+"\n"

    #Overwrite the config file
    fout = open(final_sample_name+"/"+final_sample_name+"_config.txt","w")
    for line in all_lines:
        fout.write(line)
    fout.close()
    print(all_lines)

    #Run quickvina 
    quickvina_line = '.././QuickVina2-GPU-2-1 --config '+pwd+'/'+final_sample_name+'/'+final_sample_name+'_config.txt > '+final_sample_name+'/'+final_sample_name+'_out.log\n'
    #fout = open('run_bash_script.sh','a')
    #fout.write(quickvina_line)
    #fout.close()
    os.system(quickvina_line)

    all_results = []
    drug_out_files = glob.glob(drug_folder_name+"_out/*.pdbqt")
    for drug_out_file in drug_out_files:
        fp = open(drug_out_file,"r")
        dock_score = float(fp.readlines()[1].strip("\r\n").split(":")[1].strip(" ").split(" ")[0])
        fp.close()
        drug_name = drug_out_file.split("/")[2].split("_")[0]
        all_results.append([final_sample_name,drug_name,dock_score])
    all_results_df = pd.DataFrame(all_results)
    all_results_df.columns = ['ProteinId','DrugId','Dock_Score']
    all_results_df.to_csv(final_sample_name+"/"+final_sample_name+"_dock_scores.csv",header="infer",index=None)

    #Command to clear output
    command8 = "rm -rf "+drug_folder_name+" "+drug_folder_name+"_out "
    os.system(command8)
    # -



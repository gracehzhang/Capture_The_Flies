#!/usr/bin/env python2.7
'''
    ### foldx_montecarlo_jug ###
    
    Wrapper script that runs FoldX, performing Monte 
    Carlo sampling... using jug

    Created on Jul 2, 2015

'''
__author__ = "Grace Zhang"
__version__ = "0.1"

import sys, random
import os, shutil
import argparse
from jug import TaskGenerator


@TaskGenerator
def calc_rand(obs_input,pdb_dict,filename,keep_files,i):
    trash = []
    create_accessory_files(trash,filename)
    stab_count = 0
    tot = 0
    for pdb,dict_pdb in obs_input.iteritems():
        wt = ''.join((pdb,"_repd.pdb"))
        wt_path = '/'.join((pdb_dir,wt))
        wt_copy = '/'.join((run_dir,wt))
        if os.path.isfile(wt_path):
            shutil.copyfile(wt_path, wt_copy)
            trash.append(wt_copy)
        else:
            wt0 = ''.join((pdb,".pdb"))
            wt0_path = '/'.join((pdb_dir,wt0))
            wt0_copy = '/'.join((run_dir,wt0))
            shutil.copyfile(wt0_path, wt0_copy)
            trash.append(wt0_copy)
            repairPDB(wt0)
        for m1, dict_m1 in dict_pdb.iteritems():
            m1_pdb = get_mutated_pdb_filename(m1,wt)
            m1_path = '/'.join((pdb_dir,m1_pdb))
            #m1_copy = '/'.join((run_dir,m1_pdb))
            #if os.path.isfile(m1_path):
            #   shutil.copyfile(m1_path, m1_copy)
                #trash.append(m1_copy)
            m1_nrg_path = ''.join((nrg_dir,"/",m1_pdb[:-4],".txt"))
            m1_nrg_copy = ''.join((run_dir,"/",m1_pdb[:-4],".txt"))
            #if os.path.isfile(m1_nrg_path):
            #    shutil.copyfile(m1_nrg_path, m1_nrg_copy)
                #trash.append(m1_nrg_copy)
            if not os.path.isfile(m1_path) or not os.path.isfile(m1_nrg_path):
                mutate1(m1, wt, trash, filename)
            m1_nrg = parse_nrg(m1_nrg_copy)
            for _ in range(len(dict_m1)):
                m2 = get_random_mutation(pdb,pdb_dict,m1)
                m2_pdb = get_mutated_pdb_filename(m2,m1_pdb)
                m2_nrg_path = ''.join((nrg_dir,"/",m2_pdb[:-4],".txt"))
                m2_nrg_copy = ''.join((run_dir,"/",m2_pdb[:-4],".txt"))
                if (keep_files):
                    if os.path.isfile(m2_nrg_path):
                        shutil.copyfile(m2_nrg_path, m2_nrg_copy)
                        trash.append(m2_nrg_copy)
                    else:
                        mutate2(m2, m1_pdb, trash, filename,keep_files,False)
                else:
                    mutate2(m2, m1_pdb, trash, filename,keep_files,False)
                m2_nrg = parse_nrg(m2_nrg_copy)
                if normalize(m2_nrg - m1_nrg) < 0: 
                    stab_count += 1
                tot += 1
    empty_trash("after calc_rand", trash)
    sys.stderr.write("Random simulation #" + str(i) + " complete\n")
    return stab_count / float(tot)

def repairPDB(pdb):
    '''
    create or edit runfile in folder
    call foldx with popen
        select 3 then enter repairfile
        wait for the program to finish
    return the name of the repaired pdb
    '''
    file = open("repairfile.txt", "r")
    temp = open("temp.txt", "w")
    for line in file:
        if line[0:6:] == "<PDBS>":
            temp.write(line[0:6:] + pdb + ";\n")   
        else:
            temp.write(line)
    temp.close()
    file.close()
    os.rename("temp.txt","repairfile.txt")
    os.system(foldx + "repairfile.txt" + " > /dev/null")
    os.rename(''.join((run_dir,"/RepairPDB_",pdb[:-4],".pdb")),''.join((run_dir,"/",pdb[:-4],"_repd.pdb")))
    shutil.copyfile(''.join((run_dir,"/",pdb[:-4],"_repd.pdb")),''.join((pdb_dir,"/",pdb[:-4],"_repd.pdb")))
    os.remove(''.join((run_dir,"/RepairPDB_",pdb[:-4],".fxout")))
    os.remove(''.join((run_dir,"/Unrecognized_molecules.txt")))
    
def get_random_mutation(pdb,pdb_dict,m1):
    """
    Generate a random mutation avoiding position of m1.
    """
    chain = m1[1]
    m1_pos = m1[2:-1]
    pos = [p for p in pdb_dict[pdb][chain]]
    pos.remove(m1_pos)
    rand_pos = random.choice(pos)
    
    res = [r for r in aa_table().iterkeys()]
    curr_res = pdb_dict[pdb][chain][rand_pos]
    res.remove(curr_res)
    rand_res = random.choice(res)
    return ''.join((curr_res,chain,rand_pos,rand_res))
    
def get_pdb_dict(obs_input):
    """
    Generate dictionary that contains all positions for
    each chain of each pdb file in obs_input.
    """
    pdb_dict = {}#defaultdict(lambda : defaultdict(dict))
    for pdb in obs_input:
        pdb_file = ''.join((pdb_dir,'/',pdb,'.pdb'))
        for line in open(pdb_file,'r'):
            pdb_dict.setdefault(pdb, {})
            line = line.rstrip()
            if line[0:4] == 'ATOM':
                chain = line[21]
                pos = line[22:26].lstrip()
                res = line[17:20]
                pdb_dict[pdb].setdefault(chain,{})
                pdb_dict[pdb][chain][pos] = aa_table_inv()[res]            
    return pdb_dict
      

def aa_table():
    return {"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU","G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE","P":"PRO","S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL"}
def aa_table_inv():
    return {v: k for k, v in aa_table().items()}
       
def create_accessory_files(trash, filename):
    pos_path = '/'.join((run_dir, filename))
    file = open(pos_path, "w")
    file.write("<TITLE>FOLDX_runscript;\n<JOBSTART>#;\n<PDBS>;\n<BATCH>#;\n<COMMANDS>FOLDX_commandfile;\n<PositionScan>;\n<END>#;\n<OPTIONS>FOLDX_optionfile;\n<Temperature>#;\n<pH>#;\n<END>#;\n<JOBEND>#;\n<ENDFILE>#;")
    file.close()
    trash.append(pos_path)
    rep_path = '/'.join((run_dir,"repairfile.txt"))
    file = open(rep_path, "w")
    file.write("<TITLE>FOLDX_runscript;\n<JOBSTART>#;\n<PDBS>;\n<BATCH>#;\n<COMMANDS>FOLDX_commandfile;\n<RepairPDB>#;\n<END>#;\n<OPTIONS>FOLDX_optionfile;\n<Temperature>#;\n<pH>#;\n<END>#;\n<JOBEND>#;\n<ENDFILE>#;")
    file.close()
    trash.append(rep_path)

def mutate2(m2, m1_pdb, trash, filename,keep_files,obs):
    '''
    mutate residue 2
    '''
    
    print "mutate2: " + m2
    file = open(filename, "r")
    temp = open("temp.txt", "w")
    for line in file:
        if line[0:6:] == "<PDBS>":
            temp.write(line[0:6:] + m1_pdb + ";\n") 
        elif line[0:14:] == "<PositionScan>":
            temp.write(line[0:14:] + "to_delete.txt")
            temp.write("," + m2)
            temp.write(";\n")
        else:
            temp.write(line)
    temp.close()
    file.close()
    os.rename("temp.txt", filename)
    os.system(foldx + filename + " > /dev/null")
    
    m2_mut_new_name = ''.join((m1_pdb[:-4],"_",m2,".pdb"))
    
    nrg = ''.join(("energies_",m2[2:-1],"_",m1_pdb[:-4],".txt"))
    nrg_new_name = ''.join((m2_mut_new_name[:-4],".txt"))
    nrg_path = ''.join((run_dir,'/',nrg))
    if (keep_files or obs):
        nrg_copy = ''.join((nrg_dir,'/',nrg_new_name))
        shutil.copyfile(nrg_path,nrg_copy)
    nrg_rename = ''.join((run_dir,'/',nrg_new_name)) 
    while(True):
        try:
            os.rename(nrg_path,nrg_rename)  
        except OSError:
            os.system(foldx + filename + " > /dev/null")
        else:
            break

    trash.append(nrg_rename)
    
    if os.path.isfile("to_delete.txt"):
        os.remove("to_delete.txt")
    if os.path.isfile(''.join((run_dir,"/",aa_table()[m2[0]],m2[2:-1],"_",m1_pdb))):
        os.remove(''.join((run_dir,"/",aa_table()[m2[0]],m2[2:-1],"_",m1_pdb)))
        
    m2_mut = ''.join((aa_table()[m2[-1]],m2[2:-1],"_",m1_pdb))
    
    m2_path = ''.join((run_dir,'/',m2_mut))
    if (keep_files or obs):
        m2_copy = ''.join((pdb_dir,'/',m2_mut_new_name))
        shutil.copyfile(m2_path,m2_copy)
    m2_rename = ''.join((run_dir,'/',m2_mut_new_name))
    try:
        os.rename(m2_path,m2_rename)
    except OSError:
        sys.stderr(m2_path + " not Found")
    trash.append(m2_rename)
    
   
def parse_nrg(file):
    """
    Parse deltaG value from FoldX energies file. 
    """
    nrg = 0.0
    for line in open(file,'r'):
        split = line.rstrip().split('\t')
        nrg = float(split[1])
    return nrg

def normalize(nrg):
    """
    Normalize energy value: nrg - 0.39 / 0.57
    """
    return (nrg - .39) / .57
        

def mutate1(m1, wt, trash, filename):
    '''
    create or edit positionscan file in folder
        input the correct repaired pdb
        input the correct mutation desired
    call foldx
        wait for the program to finish
    find the correct energies folder and extract the 
    '''
    file = open(filename, "r")
    temp = open("temp.txt", "w")
    for line in file:
        if line[0:6:] == "<PDBS>":
            temp.write(line[0:6:] + wt + ";\n") 
        elif line[0:14:] == "<PositionScan>":
            temp.write(line[0:14:] + "to_delete.txt," + m1 + ";\n")
        else:
            temp.write(line)
    temp.close()
    file.close()
    os.rename("temp.txt",filename)
    os.system(foldx + filename + " > /dev/null")
    if os.path.isfile(''.join((run_dir,"/",aa_table()[m1[0]],m1[2:-1],"_",wt))):
        print ''.join((run_dir,"/",aa_table()[m1[0]],m1[2:-1],"_",wt)) + " found"
        os.remove(''.join((run_dir,"/",aa_table()[m1[0]],m1[2:-1],"_",wt)))
    m1_mut = ''.join((aa_table()[m1[-1]],m1[2:-1],"_",wt))
    m1_mut_new_name = ''.join((wt[:-4],"_",m1,".pdb"))
    m1_path = ''.join((run_dir,'/',m1_mut))
    m1_copy = ''.join((pdb_dir,'/',m1_mut_new_name))
    shutil.copyfile(m1_path,m1_copy)
    m1_rename = ''.join((run_dir,'/',m1_mut_new_name))
    os.rename(m1_path,m1_rename)
    #trash.append(m1_rename)
    nrg = ''.join(("energies_",m1[2:-1],"_",wt[:-4],".txt"))
    nrg_new_name = ''.join((m1_mut_new_name[:-4],".txt"))
    nrg_path = ''.join((run_dir,'/',nrg))
    nrg_copy = ''.join((nrg_dir,'/',nrg_new_name))
    shutil.copyfile(nrg_path,nrg_copy)
    nrg_rename = ''.join((run_dir,'/',nrg_new_name))
    os.rename(nrg_path,nrg_rename)
    #trash.append(nrg_rename)

def get_mutated_pdb_filename(m, wt):
    """
    Return PDB mutated file name.
    """
    return wt[:-4] + "_" + m + ".pdb"
    
def parse_input(file):
    """
    Parse input file containing observed mutation pairs
    and store in a dictionary. The input file should be
    a tab-delimited file, with the PDB file ID in the 
    first column, mutation 1 in the second column, and 
    mutation 2 in the third column.
    """
    input_dict = {}#defaultdict(lambda : defaultdict(dict)) # Initialize dictionary of dictionaries of dictionaries
    for line in file: # For each line in file
        splt = line.rstrip().split('\t') # remove white spaces on the right side, then split at tabs.
        input_dict.setdefault(splt[0],{})
        input_dict[splt[0]].setdefault(splt[1],{})
        input_dict[splt[0]][splt[1]][splt[2]] = 1
        
    return input_dict

def empty_trash(loc, trash):
    sys.stderr.write("Cleaning up {}\n".format(loc))
    #global trash
    for t in trash:
        if os.path.isfile(t):
            os.remove(t) # empty trash
        

def parse_args(args):
    """
    Sub-function that uses argparse to parse command-line arguments.
    """
    
    parser = argparse.ArgumentParser(description=__doc__,
                                    epilog='foldx_montecarlo.py v{} written by {}'.format(__version__,__author__))
    parser.add_argument('input',
        type=argparse.FileType('r'),
        default=sys.stdin,
        help='Input file containing observed mutation')

    parser.add_argument('-n','--nsim',
        type=int,
        default=100,
        help='Number of simulations.')
    
    parser.add_argument('-keep_file',
        type=bool,
        default=False,
        help='Keep m2 files?')
    
    parser.add_argument('-obs',
        type=int,
        default=0.764705882353,
        help='Prop. in obs_input')
    
    pargs = parser.parse_args()
    
    return pargs


######## Set paths ########
foldx = "./foldx3b6 -runfile " #figure out usr/local later
run_dir = "/Users/linzhang/Desktop/CS1B class/Monte_Carlo_Test_v2/default"
pdb_dir = "/Users/linzhang/Desktop/CS1B class/Monte_Carlo_Test_v2/default/pdb"
nrg_dir = "/Users/linzhang/Desktop/CS1B class/Monte_Carlo_Test_v2/default/nrg"

###########################

######## Global variables ########
#trash = []
rand = []
##################################

if len(sys.argv) == 1:
    print "For help, run foldx_montecarlo.py with -h flag."
    sys.exit(1)

########## Get command-line arguments ##########
pargs = parse_args(sys.argv)
################################################



obs_input = parse_input(pargs.input)
pdb_dict = get_pdb_dict(obs_input)

obs = pargs.obs



'''
sys.stderr.write("error message\n")
for error log
'''

for i in range(pargs.nsim):
    filename = "positionscan" + str(i) + ".txt"
    rand_sim = calc_rand(obs_input,pdb_dict,filename,pargs.keep_file,i)
    rand.append(rand_sim)
    
   
    

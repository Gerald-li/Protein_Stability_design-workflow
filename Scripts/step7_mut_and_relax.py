#coding:utf-8
import os
import sys
#import matplotlib.pyplot as plt
import numpy as np
from pymol import cmd
from pymol import stored

def checkinputpdb(mutepdb, mutefile, relaxrange=8):
    if sys.argv[1] == "-h":
        print("USAGE: python writemutantrelax.py mutepdb mutefile relaxrange cpunumber")
    elif (os.path.isfile(mutepdb) != True):
        print("Please put original pre-mutepdb into the folder.")
        os._exit(0)
    elif (os.path.isfile(mutefile) != True):
        print("Please put mutefile into the folder.")
        os._exit(0)
    elif (relaxrange == False):
        print(relaxrange)
        print("Please input relaxrange.")
        os._exit(0)
    else:
        return(1)

def dir_exists(dir):
    if os.path.exists(dir):
        print(dir)
    if not os.path.exists(dir):
        os.mkdir(dir)

def getbestscorepdb(scorepath): 
    with open(scorepath, 'r') as sc:
        i = 0
        refscore = 0
        for eachline in sc:
            print(eachline)
            i = i + 1
            if i > 2:
                eachdata = eachline.split()
                tmpscore = float(eachdata[1])
                if tmpscore < refscore:
                    refscore = tmpscore
                    bestpdb = eachdata[21]
    return (refscore, bestpdb)

def get_bestscore_pdb_file(prepared_outfile):
    with open(prepared_outfile, 'r') as f:
        lines = f.readlines()[0].strip()
    f.close()
    bestscore_pdbname = lines.split(' ')[-1]
    bestscore_wt_relax = lines.split(' ')[0]
    return bestscore_pdbname, bestscore_wt_relax

def writemutesh(scoredpdb, mutefile):
    with open(mutefile, 'r') as mutf:
        with open("mute.sh", 'w') as mutsh:
            for eachline in mutf:
                eachline=eachline.strip("\n")
                mutsh.write("mutate.default.linuxgccrelease -s " + scoredpdb + " -mutate:mutation " + eachline + "\n")
                ##for example: mutate.default.linuxgccrelease -s the_best_score_pdb.pdb -mutate:mutation V100A
        mutsh.close()
    mutf.close()
    outfile = scoredpdb.split('.pdb')[0] + '_' + eachline + '.pdb'
    return outfile

def get_mutepdb(scoredpdb, mutation):
    with open("mute.sh", 'w') as mutsh:
        mutsh.write("mutate.default.linuxgccrelease -s ../refine/" + scoredpdb + " -mutate:mutation " + mutation + "\n")
        ##for example: mutate.default.linuxgccrelease -s the_best_score_pdb.pdb -mutate:mutation V100A
    mutsh.close()
    outfile = scoredpdb.split('.pdb')[0] + '_' + mutation + '.pdb'
    os.system('sh mute.sh')
    return outfile

def refineflag(nstruct, relax_script, outpath):  
    with open('refineflag', 'w') as f:
        f.write('-nstruct ' + str(nstruct) + '\n')   
        f.write('-relax:' + relax_script + '\n')     
        f.write('-out:path:pdb ' + outpath + '\n')   
        f.write('-out:path:score ' + outpath + '\n') 
    f.close()
    return 'refineflag'

def writerelaxsh(scoredpdb, bestscore_wt_relax, mutefile, refineflag_parameters, relaxrange, mpi_run, columns):
    with open(mutefile, 'r') as mutf:
        with open('./selected_bestscorepdb/bestscoreandpdb.csv', 'a') as fb:
            fb.write("WT_relax\t%s\t%s\n" %(bestscore_wt_relax, scoredpdb.split('.pdb')[0]))
        fb.close()
        # with open("relax.sh", 'w') as relaxsh:
        for eachline in mutf:
            eachline=eachline.strip("\n")
            print(eachline)
            dir_exists(eachline)
            os.chdir(eachline)
            refineflag(refineflag_parameters[0], refineflag_parameters[1], refineflag_parameters[2])
            dir_exists('refine')
            mutepdb = get_mutepdb(scoredpdb, eachline)
            if relaxrange != 'all':
                cmd.load(mutepdb)
                #cmd.load(mutepdb.strip(".pdb") + "_" + eachline + ".pdb", "mutepdb")
                resid=eachline[1:-1]
                # resid=resid.strip("P")
                #print(resid)
                movemap_file = eachline + "_movemap.f"
                with open(movemap_file, 'w') as mmp:
                    cmd.select("selepro", "resid " + resid)
                    cmd.select("seleproexpand", "(byres (selepro expand " + relaxrange + "))", enable=1)
                    stored.residandCAatom = []
                    cmd.iterate("seleproexpand and name CA", "stored.residandCAatom.append((resv, ID, resn))")
                    for eachdata in stored.residandCAatom:
                        print(eachdata[0])
                        mmp.write("RESIDUE " + str(eachdata[0]) + " BBCHI\n")
                mmp.close()
                # with open("relax.sh", 'w') as relaxsh:
                #     relaxsh.write("mpirun -np " + mpi_run + " relax.mpi.linuxgccrelease -s " + mutepdb.strip(
                #         ".pdb") + "_" + eachline + ".pdb -out:path:pdb ./" + eachline + " -out:path:score ./" + eachline + " @refineflag -relax:script relaxscript -in:file:movemap " + movemap_file + "\n")
                #     relaxsh.close()
                with open("relax.sh", 'w') as relaxsh:
                    relaxsh.write("mpirun --hostfile /usr/local/bin/hostfile -np " + mpi_run + " --allow-run-as-root relax.mpi.linuxgccrelease -s " + mutepdb + " @refineflag\n")
                relaxsh.close()
            if relaxrange == 'all':
                with open("relax.sh", 'w') as relaxsh:
                    relaxsh.write("mpirun --hostfile /usr/local/bin/hostfile -np " + mpi_run + " --allow-run-as-root relax.mpi.linuxgccrelease -s " + mutepdb + " @refineflag\n")
                relaxsh.close()
            os.system('sh relax.sh')
            score, bestpdb = getbestscorepdb('./refine/score.sc')
            with open('mutbestscore.out', 'w') as fm:
                fm.write("%s %s\n" %(score, bestpdb))
            fm.close()
            bestscore_pdbname, bestscore_relax = get_bestscore_pdb_file('mutbestscore.out')
            os.system('cp ./refine/' + bestscore_pdbname + '.pdb ' + '../selected_bestscorepdb/')
            #plt_distribution(eachline, './refine/score.sc', columns)
            #os.system('cp ' + eachline + '.png ' + '../selected_bestscorepdb/')
            with open('../selected_bestscorepdb/bestscoreandpdb.csv', 'a') as fb:
                fb.write("%s\t%s\t%s\n" %(eachline, score, bestpdb))
            os.chdir('../')
        fb.close()
    mutf.close()


if __name__ == '__main__':
    basedir = os.path.join(os.path.dirname(os.path.abspath(__file__)), './')
    os.chdir(basedir)
    pdbfilename = sys.argv[1]
    nstruct = sys.argv[2]
    relaxrange = sys.argv[3]
    cpunumber = sys.argv[4]
    workdir = pdbfilename.split('.pdb')[0]
    dir_exists(workdir)
    os.chdir(workdir)
    dir_exists('selected_bestscorepdb')
    prepared_outfile = 'prepared.out'
    bestscore_pdbname, bestscore_wt_relax = get_bestscore_pdb_file(prepared_outfile)
    refined_native_pdb = bestscore_pdbname + '.pdb'
    os.system('cp ./refine/%s selected_bestscorepdb' % (refined_native_pdb))
    mute_list_file ='../mut_for_relax_file.csv'
    refineflag_parameters = [nstruct, "default_repeats 5", "./refine"] 
    writerelaxsh(refined_native_pdb, bestscore_wt_relax, mute_list_file, refineflag_parameters, relaxrange, mpi_run=cpunumber, columns=9)

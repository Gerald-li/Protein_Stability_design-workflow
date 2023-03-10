#coding:utf-8
import os
import sys
#import matplotlib.pyplot as plt
import numpy as np
from pymol import cmd
from pymol import stored

def checkinputpdb(mutepdb, mutefile, relaxrange=8): #读取WT的pdb，突变列表，以及突变后准备优化的以突变残基为中心的半径范围
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
    #        shutil.rmtree(dir)
    if not os.path.exists(dir):
        os.mkdir(dir)

def getbestscorepdb(scorepath): #从refine后的score中选择能量最低的
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

def refineflag(nstruct, relax_script, outpath):  ##写出refineflag文件，供rosetta relax使用
    with open('refineflag', 'w') as f:
        f.write('-nstruct ' + str(nstruct) + '\n')          # -nstruct 代表进行几次relax计算
        f.write('-relax:' + relax_script + '\n')      #-relax:default_repeats 代表relax过程中，算法进行多少次退火模拟
        f.write('-out:path:pdb ' + outpath + '\n')    #-out:path:pdb 代表在哪个文件夹中输出结果文件(pdb格式)
        f.write('-out:path:score ' + outpath + '\n')  #-out:path:score 代表在哪个文件夹中输出打分结果文件
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

# def plt_distribution(figname, scorefile, columns): ##打分文件和需要分成几个柱状图来显示
#     total_socre = []
#     with open(scorefile, 'r') as fs:
#         lines = fs.readlines()
#     fs.close()
#     for line in lines[2:]:
#         total_socre.append(line.split()[1])
#     total_socre_np = np.array(total_socre, dtype=float)
#     plt.hist(x=total_socre_np, bins=columns)
#     plt.savefig(figname + '.png')
#     plt.clf()



#def main():
#    mutepdb = sys.argv[1]
#    mutefile = sys.argv[2]
#    relaxrange = sys.argv[3]  ##relaxrange = 'all'表示relax整个蛋白; 或者一个范围例如8，表示突变残基8A范围内的残基进行relax.

#    if checkinputpdb(mutepdb, mutefile, relaxrange):
#        print("checked")
#        writemutesh(mutepdb, mutefile)
#        # if relaxrange == 'all':
#        #     writerelax_all_sh(mutepdb, mutefile)
#        # if relaxrange != 'all':
#        writerelaxsh(mutepdb, mutefile, relaxrange=8)


if __name__ == '__main__':
    basedir = os.path.join(os.path.dirname(os.path.abspath(__file__)), './')
    os.chdir(basedir)
    # pdbfilename = '1nf2_chainA_delMg.pdb'
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
    # relaxrange = '8' ##relaxrange = 'all'，默认全蛋白优化，如果蛋白较大，建议局部优化，不小于8A（relaxrange='8'），
    # mutepdb = writemutesh(refined_native_pdb, mute_list_file)
    # refineflag_parameters = ["nstruct=50", "relax_script='default_repeats 5'", "outpath='./refine'"] ##larger nstruct means we can get a precise result, but more time cost. 20-50 recommend.
    # nstruct = "3"
    refineflag_parameters = [nstruct, "default_repeats 5", "./refine"] ##larger nstruct means we can get a precise result, but more time cost. 20-50 recommend.
    writerelaxsh(refined_native_pdb, bestscore_wt_relax, mute_list_file, refineflag_parameters, relaxrange, mpi_run=cpunumber, columns=9) #mpi_run='39'表明用39个核来运行此任务，默认为20个核。

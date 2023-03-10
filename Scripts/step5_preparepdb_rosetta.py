#coding:utf-8
import os
import sys

#pdbfilename = sys.argv[1]
def dir_exists(dir):
    if os.path.exists(dir):
        print(dir)
    #        shutil.rmtree(dir)
    if not os.path.exists(dir):
        os.mkdir(dir)

def refineflag(nstruct, relax_script, outpath):  ##写出refineflag文件，供rosetta relax使用
    with open('refineflag', 'w') as f:
        f.write('-nstruct ' + str(nstruct) + '\n')          # -nstruct 代表进行几次relax计算
        f.write('-relax:' + relax_script + '\n')      #-relax:default_repeats 代表relax过程中，算法进行多少次退火模拟
        f.write('-out:path:pdb ' + outpath + '\n')    #-out:path:pdb 代表在哪个文件夹中输出结果文件(pdb格式)
        f.write('-out:path:score ' + outpath + '\n')  #-out:path:score 代表在哪个文件夹中输出打分结果文件
    f.close()


def preparepdb(pdbbeforerefine, mpi_run): #读取要进行预测的pdb结构并进行处理，之后对处理后的pdb进行refine
    with open('prep_relax.sh','w') as f:
        f.write("score_jd2.linuxgccrelease -in:file:s " + pdbbeforerefine + " -ignore_unrecognized_res -out:pdb\n")
        if int(mpi_run) >= 2:
            f.write("mpirun --hostfile /usr/local/bin/hostfile -np " + mpi_run + " --allow-run-as-root relax.mpi.linuxgccrelease -s " + pdbbeforerefine.strip(".pdb") + "_0001.pdb @refineflag")
        if int(mpi_run) == 1:
            f.write("relax.default.linuxgccrelease -s " + pdbbeforerefine.strip(".pdb") + "_0001.pdb @refineflag")
    f.close()
    os.system('sh prep_relax.sh')
#    os.system("score_jd2.linuxgccrelease -in:file:s " + pdbbeforerefine + " -ignore_unrecognized_res -out:pdb")
    # os.system("mpirun -np 20 relax.mpi.linuxgccrelease -s " + pdbbeforerefine.strip(
    #     ".pdb") + "_0001.pdb -out:path:pdb ./refine -out:path:score ./refine @refineflag -relax:script relaxscript \n")
#    os.system("relax.default.linuxgccrelease -s " + pdbbeforerefine.strip(".pdb") + "_0001.pdb @refineflag")


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


def main_bk():
    workdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '/data/jinlong/server_test/mypro_test/')
    os.chdir(workdir)
    for root,dirs,files in os.walk(workdir):
        for dir in dirs:
            if dir.startswith('5NLT'):
                os.chdir(dir)
                pdbfilename = dir + '.pdb'
                #pdbfilename = sys.argv[1]
                dir_exists('refine')
                # refineflag(nstruct=20, relax_script='default_repeats 5', outpath='./refine')
                # preparepdb(pdbfilename)
                scoreandpdb = getbestscorepdb("./refine/score.sc")
                with open ("prepared.out",'w') as po:
                    po.write("%s %s" %(scoreandpdb[0],scoreandpdb[1]))
                os.chdir('../')


if __name__ == '__main__':
    basedir = os.path.join(os.path.dirname(os.path.abspath(__file__)), './')
    os.chdir(basedir)
    # pdbfilename = '1nf2_chainA_delMg.pdb'
    pdbfilename = sys.argv[1]
    nstruct = int(sys.argv[2])
    cpunumber = sys.argv[3]
    workdir = pdbfilename.split('.pdb')[0]
    dir_exists(workdir)
    os.system('cp ' + pdbfilename + ' ' + workdir)
    os.chdir(workdir)
    dir_exists('refine')
    # refineflag(nstruct=50, relax_script='default_repeats 5', outpath='./refine')##nstruct = 50，越大越精确，但越耗时
    refineflag(nstruct, relax_script='default_repeats 5', outpath='./refine')
    preparepdb(pdbfilename, mpi_run=cpunumber) ##mpi_run='1'时，不使用mpi
    scoreandpdb = getbestscorepdb("./refine/score.sc")
    with open("prepared.out",'w') as po:
        po.write("%s %s" %(scoreandpdb[0],scoreandpdb[1]))
    os.chdir('../')

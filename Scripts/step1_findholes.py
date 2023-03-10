# coding=utf-8
import math
import os
import sys
import multiprocessing
import time

from pymol import cmd
from pymol import stored
from multiprocessing import Process  # 多任务并行处理
from multiprocessing import cpu_count  # 读取cpu的总线程数
from multiprocessing import Pool  # 多任务并行处理管理器

#Pymol is needed for this script
#Usage： python findholes.py yourpdb.pdb cpunumbers radiuscutoff
#you will get a pdb file named holes.pdb

def getmidpoint(atom1, atom2, cutoff):  # 计算两个坐标的中点距离，如果比2*cutoff大，就保留并计算两者的中点
    dis = cmd.distance(selection1="ID " + str(atom1), selection2="ID " + str(atom2))
    if dis >= 2 * cutoff:
        cmd.center(selection="(ID " + str(atom1) + ")+(ID " + str(atom2) + ")")
        mid = cmd.get_position()
        tmp = []
        tmp.append(mid)
        tmp.append("(ID " + str(atom1) + ")+(ID " + str(atom2) + ")")
        return (tmp)
    else:
        return (0)


def getnearatomindex(atomnum):  # 检索原子周围nA半径内的原子，获得原子对列表
    x = cmd.select("ID " + str(atomnum) + " expand 1.6")
    cmd.set("dot_solvent", 1)  # 设定面积为溶剂化表面积
    y = cmd.get_area(selection="ID " + str(atomnum))
    cmd.get_area(selection="all", load_b=1)  # 将溶剂化表面积放到b-factor上
    cmd.select(name="nosurfatoms", selection="b<0.5")  # 选择溶剂化表面积小于0.5平方A的原子，即蛋白内部的原子
    if (x <= 3) and (y <= 0.5):
        stored.expand = []
        stored.extend = []
        nearatomindex = []
        cmd.iterate_state(1, "(ID " + str(atomnum) + " expand 8) and nosurfatoms",
                          "stored.expand.append((ID))")  # 获得原子周围的原子
        cmd.iterate_state(1, "(ID " + str(atomnum) + " extend 8) and nosurfatoms",
                          "stored.extend.append((ID))")  # 获得原子周围与原子所在残基相连的原子
        nearatomindex.append(atomnum)
        nearatomindex.extend(list(set(stored.expand).difference(set(stored.extend))))  # 把和原子距离近，在自己残基上下相连的残基原子从周围原子里面去除
        return (nearatomindex)
    else:
        return (0)


def removeredundancypoint(midpointlist):  # 通过排序把得到的中点位置进行去沉冗化
    sortedlist = sorted(midpointlist)
    lastlist = [0, 0, 0]
    noredundancylist = []
    for each in sortedlist:
        if each[0][0] != lastlist:
            noredundancylist.append(each)
            lastlist = each[0][0]
    return (noredundancylist)


def scanradius(atomselection, cutoff):  # 找到空隙点的最大空隙半径，0.2A扫描递增
    cmd.center(selection=atomselection[1])
    for i in range(40):
        d = 0.2
        r = cutoff - d + (i * d)
        inrangeatomnum = cmd.select("center around " + str(r))
        if (inrangeatomnum != 0):
            if r > cutoff:
                atomselection.append(r)
                return (atomselection)
            else:
                break


def getholes(pdbfilename, begin, end, cutoff):  # 找到缺陷点
    cmd.load(pdbfilename, 'pdb')
    allmidpoint = []
    for i in range(end - begin + 1):  # 找到所有原子周围的范围中的原子
        nearatomindex = getnearatomindex(i + begin)
        if nearatomindex != 0:
            atom1 = begin + i
            for j in range(len(nearatomindex)):
                atom2 = nearatomindex[j]
                midlist = getmidpoint(atom1, atom2, cutoff)  # 得到所有的点-点距离大于2*cutoff值的中间点
                if midlist != 0:
                    allmidpoint.append(midlist)
                '''with open("midpoint", 'a') as midfile:
                    midfile.write(midstr)'''
    acceptmidpoint = removeredundancypoint(allmidpoint)  # 去除重复的中间点
    radius = []
    for each in acceptmidpoint:
        result = scanradius(each, cutoff=2.4)
        if result != None:
            radius.append(result)
    if radius != [None, None]:
        return (radius)


def distant(list1, list2):  # 给出两组坐标之间的间距
    dis = math.sqrt((list1[0] - list2[0]) ** 2 + (list1[1] - list2[1]) ** 2 + (list1[2] - list2[2]) ** 2)
    return (dis)


def selepoint(pointlist):  # 将得到的空隙点去沉冗化，把距离在半径内的点去掉。
    thinlist = pointlist
    for i in range(len(thinlist)):
        tmp = thinlist[0]
        for each in thinlist:
            dis = distant(tmp[0], each[0])
            if (dis <= tmp[2]) and (tmp[2] >= each[2]):
                thinlist.remove(each)
        thinlist.append(tmp)
    if len(pointlist) != len(thinlist):
        selepoint(thinlist)
    else:
        return (thinlist)


#############################################

def main():
    pdbfilename = sys.argv[1]
    cpunumber = int(sys.argv[2])
    cutoff=float(sys.argv[3]) #空隙半径，默认2.5A
    #print(pdbfilename, cpunumber)
    cmd.load(pdbfilename, 'pdb')
    totalatomnum = cmd.select("all")
    hundrodatoms = int(totalatomnum // 100)
    swich = 1
    if swich == 0:
        outlist = getholes(195, 200, cutoff=2.5)    #调试专用，测试时候用，可以选择几个原子试验，或者只探测特定点。
    elif swich == 1:
        pool = Pool(processes=cpunumber)
        result = []
        for x in range(hundrodatoms): #把蛋白分为一百个原子一段进行并行处理
            b = x * 100 + 1
            e = (x + 1) * 100
            ret = pool.apply_async(getholes, args=(pdbfilename, b, e, cutoff))
            result.append(ret)
        ret = pool.apply_async(getholes, args=(pdbfilename, hundrodatoms * 100 + 1, totalatomnum, cutoff)) #处理最后尾部不够100个的少数原子
        result.append(ret)
        pool.close()
        pool.join()
        outlist = []
        for each in result:
            outlist.extend(each.get())

    selectedpoint = selepoint(outlist)  # 去掉沉冗的点

    with open("holes.pdb", 'w') as pocket:  # 将找到的空隙点输出到一个名为holes.pdb的文件
        for i in range(len(selectedpoint)):
            pocket.write("ATOM %6d  P   PPP A%4d    %8.3f%8.3f%8.3f  1.00 %5.2f           P  \n" % (
                i + 1, i + 1, selectedpoint[i][0][0], selectedpoint[i][0][1], selectedpoint[i][0][2],
                selectedpoint[i][2]))


if __name__ == '__main__':
    main()

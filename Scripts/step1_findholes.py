# coding=utf-8
import math
import os
import sys
import multiprocessing
import time

from pymol import cmd
from pymol import stored
from multiprocessing import Process  
from multiprocessing import cpu_count
from multiprocessing import Pool

#Pymol is needed for this script
#Usage： python findholes.py yourpdb.pdb cpunumbers radiuscutoff
#you will get a pdb file named holes.pdb

def getmidpoint(atom1, atom2, cutoff):
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


def getnearatomindex(atomnum):  
    x = cmd.select("ID " + str(atomnum) + " expand 1.6")
    cmd.set("dot_solvent", 1)  
    y = cmd.get_area(selection="ID " + str(atomnum))
    cmd.get_area(selection="all", load_b=1)  
    cmd.select(name="nosurfatoms", selection="b<0.5") 
    if (x <= 3) and (y <= 0.5):
        stored.expand = []
        stored.extend = []
        nearatomindex = []
        cmd.iterate_state(1, "(ID " + str(atomnum) + " expand 8) and nosurfatoms",
                          "stored.expand.append((ID))") 
        cmd.iterate_state(1, "(ID " + str(atomnum) + " extend 8) and nosurfatoms",
                          "stored.extend.append((ID))")
        nearatomindex.append(atomnum)
        nearatomindex.extend(list(set(stored.expand).difference(set(stored.extend))))
        return (nearatomindex)
    else:
        return (0)


def removeredundancypoint(midpointlist):
    sortedlist = sorted(midpointlist)
    lastlist = [0, 0, 0]
    noredundancylist = []
    for each in sortedlist:
        if each[0][0] != lastlist:
            noredundancylist.append(each)
            lastlist = each[0][0]
    return (noredundancylist)


def scanradius(atomselection, cutoff):
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


def getholes(pdbfilename, begin, end, cutoff):
    cmd.load(pdbfilename, 'pdb')
    allmidpoint = []
    for i in range(end - begin + 1):
        nearatomindex = getnearatomindex(i + begin)
        if nearatomindex != 0:
            atom1 = begin + i
            for j in range(len(nearatomindex)):
                atom2 = nearatomindex[j]
                midlist = getmidpoint(atom1, atom2, cutoff)
                if midlist != 0:
                    allmidpoint.append(midlist)
                '''with open("midpoint", 'a') as midfile:
                    midfile.write(midstr)'''
    acceptmidpoint = removeredundancypoint(allmidpoint) 
    radius = []
    for each in acceptmidpoint:
        result = scanradius(each, cutoff=2.4)
        if result != None:
            radius.append(result)
    if radius != [None, None]:
        return (radius)


def distant(list1, list2): 
    dis = math.sqrt((list1[0] - list2[0]) ** 2 + (list1[1] - list2[1]) ** 2 + (list1[2] - list2[2]) ** 2)
    return (dis)


def selepoint(pointlist): 
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
    cutoff=float(sys.argv[3])
    #print(pdbfilename, cpunumber)
    cmd.load(pdbfilename, 'pdb')
    totalatomnum = cmd.select("all")
    hundrodatoms = int(totalatomnum // 100)
    swich = 1
    if swich == 0:
        outlist = getholes(195, 200, cutoff=2.5) 
    elif swich == 1:
        pool = Pool(processes=cpunumber)
        result = []
        for x in range(hundrodatoms):
            b = x * 100 + 1
            e = (x + 1) * 100
            ret = pool.apply_async(getholes, args=(pdbfilename, b, e, cutoff))
            result.append(ret)
        ret = pool.apply_async(getholes, args=(pdbfilename, hundrodatoms * 100 + 1, totalatomnum, cutoff)) 
        result.append(ret)
        pool.close()
        pool.join()
        outlist = []
        for each in result:
            outlist.extend(each.get())

    selectedpoint = selepoint(outlist)  

    with open("holes.pdb", 'w') as pocket:  
        for i in range(len(selectedpoint)):
            pocket.write("ATOM %6d  P   PPP A%4d    %8.3f%8.3f%8.3f  1.00 %5.2f           P  \n" % (
                i + 1, i + 1, selectedpoint[i][0][0], selectedpoint[i][0][1], selectedpoint[i][0][2],
                selectedpoint[i][2]))


if __name__ == '__main__':
    main()

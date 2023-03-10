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
import numpy as np


#USAGE： python detectholetype.py pdbfilename.pdb(protein pdb) holefilename.pdb(hole.pdb of protein) radius(rasidue around nA of holes)



class residueproperty:
    def resilist(self):
        residueoneletterlist = ["G", "A", "S", "C", "P", "T", "V", "D", "N", "I", "L", "M", "E", "Q", "K", "H", "F",
                                "Y", "W", "R"]
        return (residueoneletterlist)

    def resilength(self, letter=None):
        residuelengthdict = {"G": "0", "A": "1", "S": "2", "C": "2", "P": "2", "T": "3", "V": "3", "D": "3", "N": "3",
                             "I": "4", "L": "4", "M": "5", "E": "5", "Q": "5", "K": "6", "H": "6", "F": "7", "Y": "8",
                             "W": "10", "R": "8"}
        if letter != None:
            return (residuelengthdict[letter])

    def resicharge(self, letter=None):
        residuechargedict = {"G": "0", "A": "0", "S": "0", "C": "0", "P": "0", "T": "0", "V": "0", "D": "-1", "N": "0",
                             "I": "0", "L": "0", "M": "0", "E": "-1", "Q": "0", "K": "1", "H": "1", "F": "0", "Y": "0",
                             "W": "0", "R": "1"}
        if letter != None:
            return (residuechargedict[letter])

    def resipka(self, letter=None):
        residuepkadict = {"G": "0", "A": "0", "S": "0", "C": "8.33", "P": "0", "T": "0", "V": "0", "D": "3.86",
                          "N": "0", "I": "0", "L": "0", "M": "0", "E": "4.25", "Q": "0", "K": "10.79", "H": "6.04",
                          "F": "0", "Y": "10.07", "W": "0", "R": "12.48"}
        if letter != None:
            return (residuepkadict[letter])

    def resihbond(self, letter=None):
        residuehbonddict = {"G": "0", "A": "0", "S": "1", "C": "1", "P": "0", "T": "1", "V": "0", "D": "2", "N": "2",
                            "I": "0", "L": "0", "M": "0", "E": "2", "Q": "2", "K": "1", "H": "2", "F": "0", "Y": "1",
                            "W": "1", "R": "2"}
        if letter != None:
            return (residuehbonddict[letter])

    def resihydropho(self, letter=None):
        residuehydrophobicdict = {"G": "1.96", "A": "2.05", "S": "1.54", "C": "1.84", "P": "1.73", "T": "1.62",
                                  "V": "2.34", "D": "1.08", "N": "1.16", "I": "2.53", "L": "2.33", "M": "2.06",
                                  "E": "1.18", "Q": "1.11", "K": "0.7", "H": "1.4", "F": "2.41", "Y": "1.82",
                                  "W": "2.17", "R": "0"}
        if letter != None:
            return (residuehydrophobicdict[letter])


def getaroundresiduelist(pdbfilename, holesfilename, radius=4.0):
    cmd.load(pdbfilename, format='pdb', object='model1')
    cmd.load(holesfilename, format='pdb', object="holes")
    print(holesfilename)
    totalholes = cmd.select("holes")
    stored.nearID = []
    stored.neartype = []
    residuesaroundholes = []
    for eachpoint in range(totalholes):
        atomid = eachpoint + 1
        stored.nearID = []
        stored.neartype = []
        residuesaroundpoint = []
        selenumber = cmd.select("hole_" + str(atomid), "/holes///" + str(atomid) + " around 4 and /model1")
        print(selenumber)
        cmd.iterate("hole_" + str(atomid), "stored.nearID.append(resv)")
        cmd.iterate("hole_" + str(atomid), "stored.neartype.append(oneletter)")
        print(stored.nearID)
        print(stored.neartype)
        newID = list(zip(stored.nearID, stored.neartype))
        nnewID = list(set(newID))
        residuesaroundpoint.append(atomid)
        residuesaroundpoint.append(nnewID)
        print('!!!!!!!!!!!!!!!!!!!!')
        print(residuesaroundpoint)
        residuesaroundholes.append(residuesaroundpoint)
        # temp=cmd.iterate("hole_" + str(eachpoint + 1), "print(name,resn,resv)")  
        # cmd.save("hole_" + str(eachpoint + 1) + ".pdb", selection="hole_" + str(eachpoint + 1), format='pdb')
    return residuesaroundholes


def getmutlist(restype, length=9, charge=9, pka=9, hbond=9, hydrophobic=9):
    lengthout = []
    chargeout = []
    pkaout = []
    hbondout = []
    hydrophobicout=[]
    allout = []

    residuelengthdict = {"G": "0", "A": "1", "S": "2", "C": "2", "P": "2", "T": "3", "V": "3", "D": "3", "N": "3",
                         "I": "4", "L": "4", "M": "5", "E": "5", "Q": "5", "K": "6", "H": "6", "F": "7", "Y": "8",
                         "W": "10", "R": "8"}
    residuechargedict = {"G": "0", "A": "0", "S": "0", "C": "0", "P": "0", "T": "0", "V": "0", "D": "-1", "N": "0",
                         "I": "0", "L": "0", "M": "0", "E": "-1", "Q": "0", "K": "1", "H": "1", "F": "0", "Y": "0",
                         "W": "0", "R": "1"}
    residuepkadict = {"G": "0", "A": "0", "S": "0", "C": "8.33", "P": "0", "T": "0", "V": "0", "D": "3.86", "N": "0",
                      "I": "0", "L": "0", "M": "0", "E": "4.25", "Q": "0", "K": "10.79", "H": "6.04", "F": "0",
                      "Y": "10.07", "W": "0", "R": "12.48"}
    residuehbonddict = {"G": "0", "A": "0", "S": "1", "C": "1", "P": "0", "T": "1", "V": "0", "D": "2", "N": "2",
                        "I": "0", "L": "0", "M": "0", "E": "2", "Q": "2", "K": "1", "H": "2", "F": "0", "Y": "1",
                        "W": "1", "R": "2"}
    residuehydrophobicdict = {"G": "1.96", "A": "2.05", "S": "1.54", "C": "1.84", "P": "1.73", "T": "1.62", "V": "2.34",
                              "D": "1.08", "N": "1.16",
                              "I": "2.53", "L": "2.33", "M": "2.06", "E": "1.18", "Q": "1.11", "K": "0.7", "H": "1.4",
                              "F": "2.41", "Y": "1.82",
                              "W": "2.17", "R": "0"}
    respropertylist = {"length": length, "charge": charge, "pka": pka, "hbond": hbond, "hydrophobic": hydrophobic}
    for resproperty in respropertylist:
        if respropertylist[resproperty] != 9:
            propertyvar = eval(resproperty)
            propertyvardict = eval("residue" + resproperty + "dict")
            seleout = str(resproperty) + "out"
            print(propertyvardict)
            if propertyvar == 2:
                for eachres in residuelengthdict:
                    if float(propertyvardict[eachres]) > float(propertyvardict[restype]):
                        eval(seleout).append(eachres)
            if propertyvar == 1:
                for eachres in residuelengthdict:
                    if float(propertyvardict[eachres]) >= float(propertyvardict[restype]):
                        eval(seleout).append(eachres)
            if propertyvar == 0:
                for eachres in residuelengthdict:
                    if float(propertyvardict[eachres]) == float(propertyvardict[restype]):
                        eval(seleout).append(eachres)
            if propertyvar == -1:
                for eachres in residuelengthdict:
                    if float(propertyvardict[eachres]) <= float(propertyvardict[restype]):
                        eval(seleout).append(eachres)
            if propertyvar == -2:
                for eachres in residuelengthdict:
                    if float(propertyvardict[eachres]) < float(propertyvardict[restype]):
                        eval(seleout).append(eachres)

            print("!!!!!!!!!!!!!!!!!!!!!!")
            print(resproperty)
            print(eval(seleout))
            allout.append(eval(seleout))
    print(allout)
    intersectionout = allout[0]
    for out in allout:
        intersectionout = list(set(intersectionout).intersection(out))
    return (intersectionout)


def judgeenvioroment(pdbfilename, holesfilename, radius, ):
    rp = residueproperty()
    nearresiduelist = getaroundresiduelist(pdbfilename, holesfilename, radius)
    if os.path.exists("probe_result"):
        os.remove("probe_result")
    for eachprobe in nearresiduelist:
        lengthtmp = []
        chargetmp = []
        pkatmp = []
        hbondtmp = []
        hydrophobictmp = []
        for eachresi in eachprobe[1]:
            lengthtmp.append(rp.resilength(eachresi[1]))
            chargetmp.append(rp.resicharge(eachresi[1]))
            pkatmp.append(rp.resipka(eachresi[1]))
            hbondtmp.append(rp.resihbond(eachresi[1]))
            hydrophobictmp.append(rp.resihbond(eachresi[1]))
        print(lengthtmp)
        newlengthtmp = list(map(float, lengthtmp))
        newchargetmp = list(map(float, chargetmp))
        newpkatmp = list(map(float, pkatmp))
        newhbondtmp = list(map(float, hbondtmp))
        newhydrophobictmp = list(map(float, hydrophobictmp))

        lengthmin = np.min(newlengthtmp)
        lengthmax = np.max(newlengthtmp)
        lengthaverage = np.mean(newlengthtmp)
        chargemin = np.min(newchargetmp)
        chargemax = np.max(newchargetmp)
        chargeaverage = np.mean(newchargetmp)
        pkamin = np.min(newpkatmp)
        pkamax = np.max(newpkatmp)
        pkaaverage = np.mean(newpkatmp)
        hbondmin = np.min(newhbondtmp)
        hbondmax = np.max(newhbondtmp)
        hbondaverage = np.mean(newhbondtmp)
        hydrophobicmin = np.min(newhydrophobictmp)
        hydrophobicmax = np.max(newhydrophobictmp)
        hydrophobicaverage = np.mean(newhydrophobictmp)

        with open("probe_result", 'a') as f:
            f.write("probe number is: {:d}\n".format(eachprobe[0]))
            f.write("probe lengthproperty is min: {:f}, max: {:f}, average: {:f}\n".format(lengthmin, lengthmax,
                                                                                           lengthaverage))
            f.write("probe chargeproperty is min: {:f}, max: {:f}, average: {:f}\n".format(chargemin, chargemax,
                                                                                           chargeaverage))
            f.write("probe pkaproperty is min: {:f}, max: {:f}, average: {:f}\n".format(pkamin, pkamax,
                                                                                        pkaaverage))
            f.write("probe hbondproperty is min: {:f}, max: {:f}, average: {:f}\n".format(hbondmin, hbondmax,
                                                                                          hbondaverage))
            f.write("probe hydrophobicproperty is min: {:f}, max: {:f}, average: {:f}\n".format(hydrophobicmin,
                                                                                                hydrophobicmax,
                                                                                                hydrophobicaverage))


#    judgehydrophobicity=

# def judgeholetype():

# def givemutelist():


def main():
    pdbfilename = sys.argv[1]
    holesfilename = sys.argv[2]
    radius = float(sys.argv[3])
    probelist = getaroundresiduelist(pdbfilename, holesfilename, radius)
    if os.path.exists("mute_sugguest"):
        os.remove("mute_sugguest")
    for eachprobe in probelist:
        probenumber = eachprobe[0]
        probearoundresi = eachprobe[1]
        with open("mute_sugguest", 'a') as f:
            f.write("probe_number: {:d}\n".format(probenumber))
            for eachresi in probearoundresi:
                resiID = eachresi[0]
                resitype = eachresi[1]
                mutelist = getmutlist(resitype, length=1,
                                      charge=0, hydrophobic=2) 
                listout = ''.join(mutelist)
                print(listout)
                f.write("resiID:{:s}{:d}   mute_to:{:s}\n".format(resitype, resiID, listout))

    judgeenvioroment(pdbfilename, holesfilename, radius)


if __name__ == '__main__':
    main()

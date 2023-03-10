# coding=utf-8
import os, sys
from collections import defaultdict

def get_same_element_index(ob_list, word):
    return [i for (i,v) in enumerate(ob_list) if v == word]

def get_mut_suggest_list(CNA_local_indices_file, Rig_Ind_value, mut_suggest_file, holes_join_cna_file):
    with open(CNA_local_indices_file, 'r') as fc:
        lines_cna = fc.readlines()
        fc.close()
        fs = open('selected_residues_from_cna.csv','w')
        fs.write(lines_cna[0][10:17].strip() + ',' + lines_cna[0][62:70].strip() + '\n')
        for line in lines_cna[1:]:
            resi_cna = line[10:17].strip()
            Rig_Ind = line[62:70].strip()
            if float(Rig_Ind) > Rig_Ind_value:
                fs.write(resi_cna + ',' + Rig_Ind + '\n')
        fs.close()
    with open(mut_suggest_file, 'r') as f:
        lines_mute = f.readlines()
        f.close()
    d = defaultdict(list)
    join_list = []
    resi_list = []
    with open('selected_residues_from_cna.csv','r') as f1:
        lines = f1.readlines()
        for line in lines[1:]:
            resi = line.strip().split(',')[0]
            for mut in lines_mute:
                if mut.startswith('resiID') and not mut.strip().endswith(':'):
                    m = mut.strip().split(' ')[0].split(':')[-1]
                    if resi == m[1:]:
                        if resi in resi_list:
                            print (resi)
                        if resi not in resi_list:
                            resi_list.append(resi)
                            select_mut_index = get_same_element_index(lines_mute, mut)
                            print (select_mut_index)
                            for mut_index in select_mut_index:
                                probes=[]
                                for i in lines_mute[0:int(mut_index)+1]:
                                    if i.startswith('probe_number'):
                                        probes.append(i.strip())
                                join_list.append(mut.strip())
                                d[mut.strip()].append(probes[-1])
        # print d.values()
    with open(holes_join_cna_file,'w') as f2:
        for key, value in d.items():
            f2.write(key + ' --> ' + ','.join(value) + '\n')
        # f2.write('\n'.join(set(join_list)))
        f2.close()
    with open('holes.pdb', 'r') as fh:
        lines_holes = fh.readlines()
    fh.close()
    with open('selected_holes.pdb', 'w') as fs:
        probe_id = []
        for value in sorted(d.values()):
            # print value
            for i in value:
                j = i.split(':')[-1].strip()
                # if j not in probe_id:
                probe_id.append(int(j))
        for id in sorted(set(probe_id)):
            fs.write(lines_holes[id-1])
    with open(holes_join_cna_file, 'r') as fh:
        lines_fh = fh.readlines()
        fm = open('mut_for_relax_file.csv', 'w')
        for line in lines_fh:
            mutresidue = line.split('-->')[0].strip().split(':')[-1]
            nativeresidue = line.split('-->')[0].strip().split('mute_to:')[0].strip().split(':')[-1]
            for i in range(len(mutresidue)):
                fm.write(nativeresidue + mutresidue[i] + '\n')


if __name__ == "__main__":
    CNA_local_indices_file = sys.argv[1]
    get_mut_suggest_list(CNA_local_indices_file, Rig_Ind_value=-3, mut_suggest_file='mute_sugguest',  holes_join_cna_file='holes_join_cna.csv')

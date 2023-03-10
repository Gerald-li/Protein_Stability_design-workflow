import os, sys

def get_wt_mut(mutefile, relaxrange):
    os.system('mv ' + mutefile + ' ' + mutefile + '.backup')
    with open(mutefile+'.backup', 'r') as f:
        lines = f.readlines()
    f.close()
    WT_refine = []
    with open(mutefile, 'a') as fo:
        for eachline in lines:
            eachline=eachline.strip("\n")
            print(eachline)
            if relaxrange != 'all' and eachline[:-1] not in WT_refine:
                wt_line = eachline[:-1] + eachline[0]
                fo.write(wt_line + '\n')
                fo.write(eachline + '\n')
            if relaxrange != 'all' and eachline[:-1] in WT_refine:
                fo.write(eachline + '\n')
            WT_refine.append(eachline[:-1])
    fo.close()


if __name__=='__main__':
    mutefile = 'mut_for_relax_file.csv'
    relaxrange = sys.argv[1]
    get_wt_mut(mutefile, relaxrange)

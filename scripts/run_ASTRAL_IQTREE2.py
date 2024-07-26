import sys
import glob
import os
import multiprocessing

pre = sys.argv[1]
cores = int(sys.argv[2])

pwd = os.path.abspath('./')+'/'
OGs = glob.glob(pwd+'genesortR/iqtree2/*.fa')
models = {}

def run_iqtree2_MFP(fa):
    cmd = 'iqtree2 -m MFP -B 1000 -T 1 -s '+fa
    try:
        if os.path.getsize(fa+'.contree') < 10:
            raise Exception
    except:
        os.system(cmd)

with multiprocessing.Pool(processes=cores) as pool:
    pool.map(run_iqtree2_MFP,OGs)

with open(pre+'.genetree_IQTREE2.tre','w') as tmp_out:
    for OG in OGs:
        if os.path.isfile(OG+'.contree') is True:
            with open(OG+'.contree') as tmp_file:
                tmp_out.write(tmp_file.read().strip('\n')+'\n')
        else:
            print('Please check the gene tree of '+OG)
            sys.exit()

os.system('astral -i '+pre+'.genetree_IQTREE2.tre -t '+str(cores)+' -o '+pre+'.astral_IQTREE2')


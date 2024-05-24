import sys
import glob
import os
import multiprocessing

pre = sys.argv[1]
cores = int(sys.argv[2])

pwd = os.path.abspath('./')+'/'
OGs = glob.glob(pwd+'genesortR/iqtree2/*.fa')
models = {}
fa_models = []
with open(pwd+'/partition_data.new.txt.best_model.nex') as file:
    for line in file:
        if 'charset' not in line and 'OG' in line:
            line = line.strip(',\n')
            #Q.plant+I{0.144108}+R4{0.376699,0.180167,0.239364,0.865359,0.168895,2.13099,0.0709345,5.14675}: OG0005629{7.81994},
            name = line.split(': ')[1].split('{')[0]
            model = line.lstrip(' ').split(': ')[0]
            models[name] = model
            fa_models.append(pwd+'genesortR/iqtree2/'+name+'_pruned.fa@@'+model)

def run_iqtree2(fa_model):
    fa = fa_model.split('@@')[0]
    model = fa_model.split('@@')[1]
    cmd = 'iqtree2 -m '+model+' -B 1000 -T 1 -s '+fa
    try:
        if os.path.getsize(fa+'.contree') < 10:
            raise Exception
    except:
        os.system(cmd)

with multiprocessing.Pool(processes=cores) as pool:
    pool.map(run_iqtree2, fa_models)

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


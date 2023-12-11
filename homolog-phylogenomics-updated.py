import sys
import os
import time
import glob
import multiprocessing
from datetime import timedelta
import subprocess

#Description: This script is used to predict coding-genes based on homologs, which could be further used to construct phylogenomic trees.
#Dependencies:python, perl, cd-hit, transdecoder, java, miniprot, orthofinder, FastTree, IQTREE2, MAFFT, uniqHaplo, HmmCleaner, BMGE, AlignmentCompare 
#input_files: a protein database, DNA-level fastas, RNA-level fastas
#Usage: python3 ~/Scripts/homolog-phylogenomics.py prefix database

if len(sys.argv) < 2:
    print('#Description: This script is used to predict coding-genes based on homologs, which could be further used to construct phylogenomic trees.')
    print('#Usage: python3 ~/Scripts/homolog-phylogenomics.py prefix database (optional, required if genomes or transcripts)')
    print('        run it with a folder "raw"')
    sys.exit()

#Root settings
path = sys.argv[0].rsplit('/',1)[0] #absolute path is suggested.
uniqHaplo_path = path+'/dependencies/uniqHaplo.pl'
AlignmentCompare_path = path+'/dependencies'
BMGE_path = path+'/dependencies/BMGE-1.12/BMGE.jar'


#input_files
pwd = os.path.abspath('./')
prefix = sys.argv[1]
cores = 16 #change it accordingly
occupancy = float(2/3) #change it accordingly

#log initiation
today = time.strftime("%Y-%m-%d", time.localtime())
out_log = open(pwd+'/homolog-phylogenomics.'+prefix+'.'+str(today)+'.log','w')
now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
out_log.write(str(now)+'\n')
out_log.write('Command: python3 '+' '.join(sys.argv)+'\n')
out_log.flush()

#check dependencies

out_log.write('\n#check dependencies\n')
def check_program_exist(program):
    try:
        result = subprocess.run(['which',program],capture_output=True,text=True)
        if result.returncode == 0:
            return True
        else:
            return False
    except Exception:
        return False

run_lists = ['miniprot','gff3_file_to_proteins.pl','cd-hit','orthofinder','perl','mafft','java','phylopypruner','FastTreeMP','iqtree2']
missings = []
exists = []
for run in run_lists:
    if check_program_exist(run):
        exists.append(run)
    else:
        missings.append(run)
if len(exists) == len(run_lists):
    out_log.write('All dependencies have been checked. We are going to next step.')
else:
    out_log.write('It seems to some dependencies missed: '+', '.join(missings)+'\n')
    out_log.write('Please check it/them at first before running\n')
    sys.exit()

#number check
try:
    fastas = glob.glob(pwd+'/raw/*')
except:
    print('Please provide the absolute path of fold with genomes')
    sys.exit()

species_list = []
genomes = glob.glob(pwd+'/raw/*.genomic.fasta')
transcripts = glob.glob(pwd+'/raw/*.transcript.fasta')
proteins = glob.glob(pwd+'/raw/*.pep.fasta')
num = len(set(genomes+transcripts+proteins))
out_log.write('Occupancy: '+str(round(num*occupancy))+' in '+str(num)+' genesets.\n')


#miniprot prediction and proteins sequences extraction
if len(genomes) > 0:
    try:
        database = sys.argv[2] #absolute route
        if '/' in database:
            database_name = database.split('/')[-1]
        else:
            database_name = database
    except:
        print('Please provide the absolute path of database for alignment.')
        sys.exit()
    out_log.write('\n#miniprot prediction and proteins sequences extraction from genome\n')
    now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    out_log.write(str(now)+'\n')
    out_log.flush()
    try:
        size = os.path.getsize(pwd+'/miniprot')
    except:
        os.makedirs(pwd+'/miniprot')
    for fasta in genomes:
        start_time = time.time()
        if '/' in fasta:
            fasta_name = fasta.split('/')[-1]
        else:
            fasta_name = fasta
        species_list.append(fasta_name.split('.')[0])
        gff_file = pwd+'/miniprot/'+fasta_name+'.against.'+database_name+'.gff'
        miniprot_cmd = 'miniprot -t '+str(cores)+' --gff '+fasta+' '+database+' > '+gff_file
        try:
            if os.path.getsize(gff_file) < 100000:
                raise Exception
        except:
            os.system(miniprot_cmd)
        out =open(gff_file+'3','w')
        with open(gff_file) as file:
            lines= file.readlines()
            list = []
            count = 0
            for line in lines:
                if line.startswith('#') is False and line.startswith('[') is False:
                    n = line.strip('\n').split('\t')
                    if n[2] == 'mRNA' or n[2] == 'transcript':
                        count_1 = 0
                        count +=1
                        out.write(n[0]+'\t'+'miniprot'+'\t'+'gene\t'+n[3]+'\t'+n[4]+'\t'+n[5]+'\t'+n[6]+'\t'+n[7]+'\tID=miniprot.gene.'+str(count)+'\n')
                        out.write(n[0]+'\t'+'miniprot'+'\t'+'mRNA\t'+n[3]+'\t'+n[4]+'\t'+n[5]+'\t'+n[6]+'\t'+n[7]+'\tID=miniprot.mRNA.'+str(count)+'; Parent=miniprot.gene.'+str(count)+'\n')
                    if n[2] =='CDS' or n[2] == 'exon':
                        count_1 +=1
                        out.write(n[0]+'\t'+'miniprot'+'\t'+'CDS\t'+n[3]+'\t'+n[4]+'\t'+n[5]+'\t'+n[6]+'\t'+n[7]+'\tID=miniprot.mRNA.'+str(count)+'.cds.'+str(count_1)+'; Parent=miniprot.mRNA.'+str(count)+'\n')
                        out.write(n[0]+'\t'+'miniprot'+'\t'+'exon\t'+n[3]+'\t'+n[4]+'\t'+n[5]+'\t'+n[6]+'\t'+n[7]+'\tID=miniprot.mRNA.'+str(count)+'.exon.'+str(count_1)+'; Parent=miniprot.mRNA.'+str(count)+'\n')
        out.close()

        #get proteins and remove proteins with stop '*'
        pep_file = pwd+'/miniprot/'+fasta_name.rsplit('.',1)[0]+'.miniprot.pep.fasta'
        try:
            if os.path.getsize(pep_file) < 100000:
                raise Exception
        except:
            extract_cmd = 'gff3_file_to_proteins.pl --fasta '+fasta+' --gff3 '+gff_file+'3 > '+pep_file
            os.system(extract_cmd)
        with open(pep_file) as tmp_file:
            n = tmp_file.read().split('\n>')
            with open(pwd+'/miniprot/'+fasta_name.rsplit('.',1)[0]+'.miniprot.pep.filtered.fasta','w') as tmp_out:
                for i in n:
                    i = i.lstrip('>')
                    name = i.split('\n',1)[0].split('\t')[0]
                    seq = i.split('\n',1)[1].replace('\n','')
                    if seq[:-1].count('*') < 1:
                        tmp_out.write('>'+name+'\n'+seq+'\n')
        end_time = time.time()
        spend = end_time - start_time
        out_log.write('       '+fasta_name+': '+str(timedelta(seconds=spend))+'\n')
        out_log.flush()

#TransDecoder prediction from transcripts
if len(transcripts) > 0:
    out_log.write('\n#TransDecoder prediction from transcripts\n')
    now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    out_log.write(str(now)+'\n')
    out_log.flush()
    try:
        size = os.path.getsize(pwd+'/transdecoder')
    except:
        os.makedirs(pwd+'/transdecoder')
    os.chdir(pwd+'/transdecoder')
    for fasta in transcripts:
        start_time = time.time()
        if '/' in fasta:
            fasta_name = fasta.split('/')[-1]
        else:
            fasta_name = fasta
        species_list.append(fasta_name.split('.')[0])
        try:
            size = os.path.getsize(pwd+'/transdecoder/'+fasta_name+'.transdecoder.filtered.fasta')
        except:
            try:
                size = os.path.getsize(fasta.replace('raw','transdecoder')+'.transdecoder.pep')
                os.system('ln -s '+fasta.replace('raw','transdecoder')+'.transdecoder.pep '+pwd+'/transdecoder/'+fasta_name+'.transdecoder.filtered.fasta')
            except:
                ln_cmd = 'ln -s '+fasta+' '+fasta.replace('raw','transdecoder')
                os.system(ln_cmd)
                transdecoder_cmd_1 = 'TransDecoder.LongOrfs -t '+fasta.replace('raw','transdecoder')
                transdecoder_cmd_2 = 'TransDecoder.Predict -t '+fasta.replace('raw','transdecoder')
                os.system(transdecoder_cmd_1)
                os.system(transdecoder_cmd_2)
                os.system('ln -s '+fasta.replace('raw','transdecoder')+'.transdecoder.pep '+pwd+'/transdecoder/'+fasta_name+'.transdecoder.filtered.fasta')
        end_time = time.time()
        spend = end_time - start_time
        out_log.write('       '+fasta_name+': '+str(timedelta(seconds=spend))+'\n')
        out_log.flush()

if len(proteins) > 0:
    for fasta in proteins:
        if '/' in fasta:
            fasta_name = fasta.split('/')[-1]
        else:
            fasta_name = fasta
        species_list.append(fasta_name.split('.')[0])

#cd-hit processing
out_log.write('\n#cd-hit processing\n')
now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
out_log.write(str(now)+'\n')
try:
    size = os.path.getsize(pwd+'/cd-hit')
except:
    os.makedirs(pwd+'/cd-hit')
os.chdir(pwd+'/cd-hit')
filtered_fastas = glob.glob(pwd+'/miniprot/*.filtered.fasta')+glob.glob(pwd+'/transdecoder/*.filtered.fasta')
for filtered_fasta in filtered_fastas:
    start_time = time.time()
    cd_hit_fasta = pwd+'/cd-hit/'+filtered_fasta.split('/')[-1].rsplit('.',1)[0]+'.cd-hit-0.85.fasta'
    try:
        if os.path.getsize(cd_hit_fasta) < 100000:
            raise Exception
    except:
        cd_hit_cmd = 'cd-hit -c 0.85 -T 40 -M 50000 -i '+filtered_fasta+' -o '+cd_hit_fasta
        os.system(cd_hit_cmd)
    end_time = time.time()
    spend = end_time - start_time
    out_log.write('       '+filtered_fasta.rsplit('/')[-1]+': '+str(timedelta(seconds=spend))+'\n')
    out_log.flush()


#Orthofinder processing
out_log.write('\n#Orthofinder processing\n')
now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
start_time = time.time()
out_log.write(str(now)+'\n')
out_log.flush()
try:
    size = os.path.getsize(pwd+'/'+prefix+'.orthofinder')
    os.system('rm -r '+pwd+'/'+prefix+'.orthofinder')
    os.makedirs(pwd+'/'+prefix+'.orthofinder')
except:
    os.makedirs(pwd+'/'+prefix+'.orthofinder')

os.chdir(pwd+'/'+prefix+'.orthofinder')
cd_hit_fastas = glob.glob(pwd+'/cd-hit/*.cd-hit-0.85.fasta')+proteins

full_abbr = {}
if True:
    with open('Fullname_abbr.txt','w') as out_tmp:
        list = []
        li = [chr(i) for i in range(ord('A'),ord('Z')+1)]
        tmp_list = []
        for i in li:
            tmp = i
            tmp_list.append(i)
            for n in li:
                if n not in tmp_list:
                    tmp_n = tmp+n
                    list.append(tmp_n)
        num = 0
        for cd_hit_fasta in cd_hit_fastas:
            if cd_hit_fasta.rsplit('/')[-1].split('.')[0] in set(species_list):
                input_name = cd_hit_fasta.rsplit('/')[-1]
                output_1 = input_name.rsplit('.',1)[0]+'.formated.txt'
                output_2 = input_name.rsplit('.',1)[0]+'.formated.fa'
                out_1 = open(output_1,'w')
                out_2 = open(output_2,'w')
                out_1.write(str(input_name)+'\t'+output_2+'\n')
                species_fullname = str(input_name).split('.',1)[0]
                sp_name = 'SP'+list[num]
                out_tmp.write(str.upper(sp_name)+'\t'+species_fullname+'\n')
                full_abbr[sp_name] = species_fullname
                with open(cd_hit_fasta) as tmp_file:
                    tmp_n = tmp_file.read().split('\n>')
                    count = 1
                    for tmp_i in tmp_n:
                        tmp_i = tmp_i.lstrip('>')
                        if len(tmp_i) > 50:
                            gene_name = tmp_i.split('\n',1)[0]
                            sequence = tmp_i.split('\n',1)[1]
                            gene_name_1 = str.upper(sp_name)+'|Contig'+str(count).zfill(5)
                            out_1.write(gene_name+'\t'+gene_name_1+'\n')
                            out_2.write('>'+gene_name_1+'\n'+sequence+'\n')
                        count +=1
                out_1.close()
                out_2.close()
                num+=1

os.chdir(pwd)
orthofinder_cmd = 'orthofinder -f '+prefix+'.orthofinder -og -t '+str(cores)
try:
    num_OGs = glob.glob(pwd+'/'+prefix+'.orthofinder/OrthoFinder/Results_*/Orthogroup_Sequences/*.fa')
    if len(num_OGs) < 1000:
        raise Exception
except:
    os.system('rm -rf '+pwd+'/'+prefix+'.orthofinder/OrthoFinder/')
    os.system(orthofinder_cmd)
end_time = time.time()
spend = end_time - start_time
out_log.write('\n       Running time: '+str(timedelta(seconds=spend))+'\n')
out_log.flush()

#Phylogenomic analysis
out_log.write('\n#Phylogenomic analysis\n')
now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
out_log.write(str(now)+'\n')
out_log.flush()
os.makedirs(pwd+'/'+prefix+'.Phylogenomics')
os.chdir(pwd+'/'+prefix+'.Phylogenomics')
os.makedirs(pwd+'/'+prefix+'.Phylogenomics/rejected_OGs')
OGs_wd = glob.glob(pwd+'/'+prefix+'.orthofinder/OrthoFinder/Results_*/Orthogroup_Sequences/')[0]
OG_num = 0
while OG_num < 20000:
    tmp_num = str(OG_num).zfill(7)
    OG_num += 1
    OG_name = OGs_wd+'/OG'+tmp_num+'.fa'
    with open(OG_name) as tmp_file:
        num_OGs = len(tmp_file.read().split('\n>'))
        if num_OGs < float(num*2/3):
            break
        os.system('ln -s '+OG_name+' ./')

def remove_linebreak(fa):
    with open(fa) as tmp_file:
        n = tmp_file.read().split('\n>')
        if len(n) > 2 and os.path.getsize(fa) > 0:
            with open(fa,'w') as tmp_out:
                tmp = ''
                for i in n:
                    i = i.lstrip('>')
                    tmp_name = i.split('\n',1)[0]
                    tmp_seq = i.split('\n',1)[1].replace('\n','')
                    tmp = tmp+'>'+tmp_name+'\n'+tmp_seq+'\n'
                tmp_out.write('>'+tmp.strip('\n').lstrip('>'))
        else:
            os.system('rm '+fa)

os.makedirs(pwd+'/'+prefix+'.Phylogenomics/rejected_few_taxa_1')
os.makedirs(pwd+'/'+prefix+'.Phylogenomics/backup_preUniqHaplo')
os.makedirs(pwd+'/'+prefix+'.Phylogenomics/backup_alignments')
os.makedirs(pwd+'/'+prefix+'.Phylogenomics/backup_pre-BMGE')
os.makedirs(pwd+'/'+prefix+'.Phylogenomics/rejected_few_taxa_2')

#Check taxa number (2/3): 1st
start_time = time.time()
OGs = glob.glob(pwd+'/'+prefix+'.Phylogenomics/*.fa')
for OG in OGs:
    taxa = {}
    length = []
    with open(OG) as file:
        n = file.read().split('\n>')
        if len(n) >= float(num*2/3):
            for i in n:
                i = i.lstrip('>')
                taxa_sp = i.split('|')[0]
                seq = i.split('\n',1)[1].replace('\n','')
                if taxa_sp not in taxa.keys():
                    taxa[taxa_sp] = 1
                else:
                    taxa[taxa_sp] = taxa[taxa_sp] + 1
                length.append(len(seq))
    if len(taxa) < float(num*2/3):
        os.system('mv '+OG+' '+pwd+'/'+prefix+'.Phylogenomics/rejected_few_taxa_1')
    else:
        remove_linebreak(OG)
end_time = time.time()
spend = end_time - start_time
out_log.write('       Check occupancy ('+str(round(num*2/3))+'/'+str(num)+') 1st: '+str(timedelta(seconds=spend))+'\n')
out_log.flush()

#Remove redundant sequences using uniqHaplo
start_time = time.time()
def UniqHaplq(fa):
    os.system("perl "+uniqHaplo_path+" -a "+fa+' > '+fa+'.uniq')
    os.system('rm '+fa)
    os.system('mv '+fa+'.uniq '+fa)
os.system('cp '+pwd+'/'+prefix+'.Phylogenomics/*.fa '+pwd+'/'+prefix+'.Phylogenomics/backup_preUniqHaplo')
fas_preUniqHaplo = glob.glob(pwd+'/'+prefix+'.Phylogenomics/*.fa')
with multiprocessing.Pool(processes=cores) as pool:
    pool.map(UniqHaplq, fas_preUniqHaplo)
end_time = time.time()
spend = end_time - start_time
out_log.write('       Remove redundant sequences using uniqHaplo: '+str(timedelta(seconds=spend))+'\n')
out_log.flush()

#Mafft alignment
start_time = time.time()
def mafft_aln(fa):
    os.system("mafft --auto --localpair --quiet --maxiterate 1000 "+fa+' > '+fa+'.aln')
    os.system('mv '+fa+'.aln '+fa)
fas_mafft = glob.glob('*.fa')
with multiprocessing.Pool(processes=cores) as pool:
    pool.map(mafft_aln, fas_mafft)
end_time = time.time()
spend = end_time - start_time
out_log.write('       Mafft alignment: '+str(timedelta(seconds=spend))+'\n')
out_log.flush()

#Clean alignments with HmmCleaner
def Hmmcleaner(fa):
    os.system('HmmCleaner.pl '+fa+' --specificity')
    os.system('mv '+fa+' ./HmmCleaner_files')
    os.system('mv '+fa.replace('.fa','_hmm.log')+' ./HmmCleaner_files')
    os.system('mv '+fa.replace('.fa','_hmm.score')+' ./HmmCleaner_files')
    os.system('mv '+fa.replace('.fa','_hmm.fasta')+' '+fa)
if check_program_exist('HmmCleaner.pl'):
    start_time = time.time()
    fas_hmmcleaner = glob.glob('*.fa')
    with multiprocessing.Pool(processes=cores) as pool:
        pool.map(Hmmcleaner, fas_hmmcleaner)
    end_time = time.time()
    spend = end_time - start_time
    out_log.write('       Clean alignments with HmmCleaner: '+str(timedelta(seconds=spend))+'\n')
else:
    out_log.write('       Clean alignments with HmmCleaner: HmmCleaner Not Found and skipped\n')
out_log.flush()

#Trim alignments with BMGE
start_time = time.time()
def BMGE(fa):
    os.system('java -jar '+BMGE_path+' -t AA -i '+fa+' -of '+fa+'.BMGE')
    os.system('mv '+fa+'.BMGE '+fa)
fas_BMGE = glob.glob('*.fa')
with multiprocessing.Pool(processes=cores) as pool:
    pool.map(BMGE, fas_BMGE)
end_time = time.time()
spend = end_time - start_time
out_log.write('       Trim alignments with BMGE: '+str(timedelta(seconds=spend))+'\n')
out_log.flush()

#Remove any sequences that don't overlap with all other sequences by at least 20 amino acids. This check runs until all sequences overlap with all other sequences by at least 20 amino acids.
start_time = time.time()
def AlignmentCompare(fa):
    os.system('java -cp '+AlignmentCompare_path+' AlignmentCompare '+fa)
fas_aligncompare = glob.glob('*.fa')
with multiprocessing.Pool(processes=10) as pool:
    pool.map(AlignmentCompare, fas_aligncompare)
end_time = time.time()
spend = end_time - start_time
out_log.write('       AlignmentCompare (<20 amino acids): '+str(timedelta(seconds=spend))+'\n')
out_log.flush()

#Check taxa number (2/3): 2nd
start_time = time.time()
OGs = glob.glob(pwd+'/'+prefix+'.Phylogenomics/*.fa')
for OG in OGs:
    taxa = {}
    length = []
    if os.path.getsize(OG) < 100:
        os.system('mv '+OG+' '+pwd+'/'+prefix+'.Phylogenomics/rejected_few_taxa_2')
    else:
        try:
            with open(OG) as file:
                n = file.read().split('\n>')
                if len(n) >= float(num*2/3):
                    for i in n:
                        i = i.lstrip('>')
                        taxa_sp = i.split('|')[0]
                        seq = i.split('\n',1)[1].replace('\n','')
                        if taxa_sp not in taxa.keys():
                            taxa[taxa_sp] = 1
                        else:
                            taxa[taxa_sp] = taxa[taxa_sp] + 1
                        length.append(len(seq))
            if len(taxa) < float(num*2/3) or len(seq) < 100:
                os.system('mv '+OG+' '+pwd+'/'+prefix+'.Phylogenomics/rejected_few_taxa_2')
            else:
                remove_linebreak(OG)
        except:
            #print(OG)
            os.system('mv '+OG+' '+pwd+'/'+prefix+'.Phylogenomics/rejected_few_taxa_2')
            continue
end_time = time.time()
spend = end_time - start_time
out_log.write('       Check occupancy ('+str(round(num*2/3))+'/'+str(num)+') 2nd: '+str(timedelta(seconds=spend))+'\n')
out_log.flush()

#FastTree for each OG
start_time = time.time()
fas_fasttree = glob.glob(pwd+'/'+prefix+'.Phylogenomics/*.fa')
out_fasttree = open('run_fasttree.sh','w')
out_fasttree.write('export OMP_NUM_THREADS='+str(cores)+'\n')
for fa_fasttree in fas_fasttree:
    out_fasttree.write('FastTreeMP -slow -gamma '+ fa_fasttree +' > '+fa_fasttree+'.tre\n')
    tre_name = fa_fasttree.rsplit('.fa',1)[0]+'.tre'
    out_fasttree.write('mv '+fa_fasttree+'.tre '+tre_name+'\n')
out_fasttree.close()
os.system('sh run_fasttree.sh')
end_time = time.time()
spend = end_time - start_time
out_log.write('       FastTree: '+str(timedelta(seconds=spend))+'\n')
out_log.flush()

#check the '|' in the name and remove OGs with more than two lengths
fas = glob.glob(pwd+'/'+prefix+'.Phylogenomics/*.fa')
for fa in fas:
    length = []
    with open(fa) as tmp_file:
        tmp_n = tmp_file.read().split('\n>')
        if len(tmp_n) > round(num*2/3)-1:
            for tmp_i in tmp_n:
                if '|' not in tmp_i:
                    os.system('rm '+fa)
                    os.system('rm '+fa.rsplit('.',1)[0]+'.tre')
                else:
                    length.append(len(tmp_i.split('\n',1)[1].replace('\n','')))
        if len(set(length)) !=1:
            os.system('rm '+fa)
            os.system('rm '+fa.rsplit('.',1)[0]+'.tre')

#Phylopyprunner processing
start_time = time.time()
phylopypruner_cmd = 'phylopypruner --threads '+str(cores)+' --min-taxa '+ str(round(num*2/3))+' --min-len 100 --dir . --min-support 0.75 --mask pdist --trim-divergent 0.75 --min-pdist 0.01 --prune MI'
#out_log.write('       phylopypruner: '+phylopypruner_cmd+'\n')
out_log.flush()
os.system(phylopypruner_cmd)
end_time = time.time()
spend = end_time - start_time
out_log.write('       phylopypruner: '+str(timedelta(seconds=spend))+'\n')
out_log.flush()

#FastTree for supermatrix.fas
start_time = time.time()
os.chdir(pwd+'/'+prefix+'.Phylogenomics/phylopypruner_output')
fasttree_cmd = 'FastTreeMP -slow -gamma supermatrix.fas > FastTree.tre'
#out_log.write('       FastTree for supermatrix.fas: '+fasttree_cmd+'\n')
out_log.flush()
try:
    if os.path.getsize(pwd+'/'+prefix+'.Phylogenomics/phylopypruner_output/FastTree.tre') < 10:
        raise Exception
except:
    os.system(fasttree_cmd)
with open('FastTree.tre') as tmp_file:
    fasttree_content = tmp_file.read()
    with open(pwd+'/'+prefix+'.Phylogenomics/phylopypruner_output/FastTree.full.tre','w') as fasttree_out:
        for abbr in full_abbr.keys():
            fasttree_content = fasttree_content.replace(abbr,full_abbr[abbr])
        fasttree_out.write(fasttree_content.replace(abbr,full_abbr[abbr]))
os.system('ln -s '+pwd+'/'+prefix+'.Phylogenomics/phylopypruner_output/FastTree.full.tre '+pwd+'/'+prefix+'.FastTree.full.tre')
end_time = time.time()
spend = end_time - start_time
out_log.write('       FastTree for supermatrix.fas: '+str(timedelta(seconds=spend))+'\n')
out_log.flush()

#IQ-TREE2 -m MFP construction
start_time = time.time()
iqtree2_cmd = 'iqtree2 --threads-max '+str(cores)+' -T AUTO -B 1000 -s supermatrix.fas -Q partition_data.txt -m MFP'
#out_log.write('       IQ-TREE2 -m MFP construction: '+iqtree2_cmd+'\n')
out_log.flush()
try:
    if os.path.getsize(pwd+'/'+prefix+'.Phylogenomics/phylopypruner_output/partition_data.txt.contree') < 10:
        raise Exception
except:
    os.system(iqtree2_cmd)
with open('partition_data.txt.contree') as tmp_file:
    iqtree_content = tmp_file.read()
    with open('partition_data.txt.full.contree','w') as iqtree_out:
        for abbr in full_abbr.keys():
            iqtree_content = iqtree_content.replace(abbr,full_abbr[abbr])
        iqtree_out.write(iqtree_content.replace('_',' '))
os.system('ln -s '+pwd+'/'+prefix+'.Phylogenomics/phylopypruner_output/partition_data.txt.full.contree '+pwd+'/'+prefix+'.IQTREE2.full.tre')
end_time = time.time()
spend = end_time - start_time
out_log.write('       IQ-TREE2 -m MFP construction: '+str(timedelta(seconds=spend))+'\n')
out_log.flush()

now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
out_log.write('\nAll done: '+str(now)+'\n')
out_log.write('Please check the constructed trees: '+prefix+'.FastTree.full.tre and '+prefix+'.IQTREE2.full.tre')
out_log.close()

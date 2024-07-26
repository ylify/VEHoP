#!/usr/bin/env python3

import sys
import argparse
import os
import time
import glob
import multiprocessing
import random
from datetime import timedelta
import subprocess

st = time.time()

#Check dependencies
def check_program_exist(program):
    try:
        result = subprocess.run(['which',program],capture_output=True,text=True)
        if result.returncode == 0:
            return True
        else:
            return False
    except Exception:
        return False
run_lists = ['miniprot','bmge','gff3_file_to_proteins.pl','cd-hit','orthofinder','perl','mafft','java','phylopypruner','FastTreeMP','iqtree2','trimal']
missings = []
exists = []
for run in run_lists:
    if check_program_exist(run):
        exists.append(run)
    else:
        missings.append(run)
if len(exists) == len(run_lists):
    if '-h' not in sys.argv:
        print('All dependencies have been checked. We are going to next step.')
else:
    if len(missings) ==1 and 'HmmCleaner.pl' in missings:
        print('It seems to some dependencies missed: '+', '.join(missings))
        print('HmmCleaner.pl is optional. Pipeline will start right now.')
    else:
        print('It seems to some dependencies missed: '+', '.join(missings))
        print('Please check it/them at first before running.')
        sys.exit()

if len(sys.argv) < 2:
    os.system('python3 '+' '.join(sys.argv)+' -h')

parser = argparse.ArgumentParser(description='Description: A pipeline to construct a maximum likelihood tree (based amino acids) from genomic sequences, transcripts, and proteins.')
parser.add_argument('-p', '--prefix', type=str, help='The prefix used in the output (Required)')
parser.add_argument('-t', '--threads', type=int, help='Threads used in running (default: 40).')
parser.add_argument('-i', '--input', type=str, help='Folder containing sequences for tree construction (Required, must be in the working directory, default: raw).')
parser.add_argument('-min', '--min_taxa', type=int, help='The taxon threshold in parition (default: 2/3 of the total inputs).')
parser.add_argument('-l', '--length_cutoff', type=str, help='The length threshold in parition (default: 100).')
parser.add_argument('-mode', '--mode', type=str, help='The software (FastTree or IQTREE2) to infer the gene tree in each OG; FastTree is time-efficient and IQTREE2 is more precise (default: FastTree).')
parser.add_argument('-g', '--genetic_code', type=str, help='Genetic code for proteins prediction from transcripts, which might be different with phylum, please check by "TransDecoder.LongOrfs -h" (If the parameter is given, it will adopt TransDecoder to predict coding potential in transcripts. Optional if only proteins and genomic sequences as inputs; Required if transcripts existed in inputs, default: Universal)')
parser.add_argument('-d', '--database', type=str, help='Proteins sequences for homolog prediction from genomic and transcriptional sequences, it is suggested as proteins from its/their close relatives (three organisms from the same genus, family, order, class, or phylum are suggested, from public data) (Optional if proteins as inputs; Required if genomic or transcriptional sequences existed in inputs; It must be provided with the absolute path).')

args = parser.parse_args()

#Root settings
path = sys.argv[0].rsplit('/',1)[0] #absolute path is suggested.
uniqHaplo_path = path+'/dependencies/uniqHaplo.pl'
AlignmentCompare_path = path+'/dependencies'
BMGE_path = path+'/dependencies/BMGE-1.12/BMGE.jar'

pwd = os.path.abspath('./')

try:
    if args.input:
        input_folder = args.input
        if len(glob.glob(pwd+'/'+input_folder+'/*.fasta')) < 5:
            raise Exception
    else:
        if len(glob.glob(pwd+'/raw/*')) > 5:
            input_folder = 'raw'
except:
    print('Please provide the absolute path of fold with inputs')
    sys.exit()

if args.prefix:
    prefix = args.prefix
else:
    sys.exit()

if args.threads:
    cores = args.threads
else:
    cores = 40 #change it accordingly

if args.database:
    database = args.database
    if os.path.isfile(database) is False:
        print('Please list a vaild database with absolute path.')
        sys.exit()

if args.length_cutoff:
    seq_thres = args.length_cutoff
else:
    seq_thres = 100 #change it accordingly

if args.genetic_code:
    genetic_code = args.genetic_code
    transdecoder = 'True'
else:
    transdecoder = 'False'
    genetic_code = 'Universal'

if args.mode:
    MODE = args.mode
else:
    MODE = 'FastTree'
    
##
species_list = []
genomes = glob.glob(pwd+'/'+input_folder+'/*.genomic.fasta')
transcripts = glob.glob(pwd+'/'+input_folder+'/*.transcript.fasta')
proteins = glob.glob(pwd+'/'+input_folder+'/*.pep.fasta')
num = len(set(genomes+transcripts+proteins))
if args.min_taxa:
    occupancy = float(args.min_taxa/num)
else:
    occupancy = float(2/3) #change it accordingly
wd = pwd+'/'+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.Phylogenomics/'
same = 'False'
try:
    if os.path.getsize(pwd+'/homolog-phylogenomics.'+prefix+'.'+str(num)+'_sets.input') > 0:
        with open(pwd+'/homolog-phylogenomics.'+prefix+'.'+str(num)+'_sets.input') as tmp_file:
            pre_set = set(tmp_file.read().split('\n'))
            if set(genomes+transcripts+proteins) == pre_set:
                same = 'True'
                print('The same input files were found. It will skip orthofinder if previous one existed.')
except:
    with open(pwd+'/homolog-phylogenomics.'+prefix+'.'+str(num)+'_sets.input','w') as tmp_file:
        tmp_file.write('\n'.join(set(genomes+transcripts+proteins)))

today = time.strftime("%Y-%m-%d", time.localtime())
out_log = open(pwd+'/homolog-phylogenomics.'+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.'+str(today)+'.log','w')
now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
out_log.write(str(now)+'\n')
out_log.write('Command: python3 '+' '.join(sys.argv)+'\n')
out_log.write('Occupancy: '+str(round(num*occupancy))+' in '+str(num)+' genesets.\n')
out_log.write('Threads: '+str(cores)+'\n')
out_log.write('Length_cutoff: '+str(seq_thres)+'\n')
out_log.flush()

out_log.write('\n#check dependencies\n')
if len(exists) == len(run_lists):
    out_log.write('All dependencies have been checked. We are going to next step.\n')
else:
    if len(missings) ==1 and 'HmmCleaner.pl' in missings:
        out_log.write('It seems to some dependencies missed: '+', '.join(missings)+'\n')
        out_log.write('HmmCleaner.pl is optional. Pipeline will start right now\n')
    else:
        out_log.write('It seems to some dependencies missed: '+', '.join(missings)+'\n')
        out_log.write('Please check it/them at first before running\n')
        sys.exit()

def check_point(program_same_occu):
    try:
        with open(program_same_occu.split('@@')[0]+'.checkpoint.ok') as file:
            tmp_lines = file.readlines()
            tmp_occu = tmp_lines[0].strip('\n').split(' ')[1]
            if program_same_occu.split('@@')[1] == 'True' and tmp_occu == program_same_occu.split('@@')[2]:
                return True
            else:
                return False
    except:
        return False

def write_check(program):
    with open(program+'.checkpoint.ok','w') as tmp_out:
        tmp_out.write('Occupancy: '+str(round(occupancy,2))+'\n')
        tmp_out.write('Number of Species: '+str(len(proteins+genomes+transcripts))+'\n')
        tmp_out.write('\t'.join(set(proteins+genomes+transcripts)))

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

if len(genomes+transcripts) > 0:
    try:
        if '/' in database:
            database_name = database.split('/')[-1]
        else:
            database_name = database
        if os.path.isfile(database) is False:
            raise Exception
    except:
        if len(genomes) > 0:
            print('Please provide the absolute path of database for alignment.')
            sys.exit()
        elif transdecoder == 'Flase' and len(transcripts) > 0:
            print('Please provide the absolute path of database for alignment.')
            sys.exit()
    out_log.write('\n#miniprot prediction and proteins sequences extraction from genome or transcripts\n')
    now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    out_log.write(str(now)+'\n')
    out_log.flush()
    try:
        size = os.path.getsize(pwd+'/miniprot')
    except:
        os.makedirs(pwd+'/miniprot')
    if len(transcripts) > 0 and transdecoder == 'False':
        fastas_list = genomes+transcripts
    else:
        fastas_list = genomes
    for fasta in fastas_list:
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
        if os.path.isfile(gff_file+'3') is False:
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
        spend = round(end_time - start_time)
        out_log.write('       '+fasta_name+': '+str(timedelta(seconds=spend))+'\n')
        out_log.flush()

if transdecoder == 'True' and len(transcripts) > 0:
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
                size = os.path.getsize(fasta.replace(input_folder,'transdecoder')+'.transdecoder.pep')
                os.system('ln -s '+fasta.replace(input_folder,'transdecoder')+'.transdecoder.pep '+pwd+'/transdecoder/'+fasta_name+'.transdecoder.filtered.fasta')
            except:
                ln_cmd = 'ln -s '+fasta+' '+fasta.replace(input_folder,'transdecoder')
                os.system(ln_cmd)
                transdecoder_cmd_1 = 'TransDecoder.LongOrfs --genetic_code '+genetic_code+ ' -t '+fasta.replace(input_folder,'transdecoder')
                transdecoder_cmd_2 = 'TransDecoder.Predict --genetic_code '+genetic_code+ ' -t '+fasta.replace(input_folder,'transdecoder')
                os.system(transdecoder_cmd_1)
                os.system(transdecoder_cmd_2)
                os.system('ln -s '+fasta.replace('raw','transdecoder')+'.transdecoder.pep '+pwd+'/transdecoder/'+fasta_name+'.transdecoder.filtered.fasta')
        end_time = time.time()
        spend = round(end_time - start_time)
        out_log.write('       '+fasta_name+': '+str(timedelta(seconds=spend))+'\n')
        out_log.flush()

if len(proteins) > 0:
    for fasta in proteins:
        if '/' in fasta:
            fasta_name = fasta.split('/')[-1]
        else:
            fasta_name = fasta
        species_list.append(fasta_name.split('.')[0])

if len(genomes+transcripts) > 0:
    out_log.write('\n#cd-hit processing\n')
    now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    out_log.write(str(now)+'\n')
    try:
        size = os.path.getsize(pwd+'/cd-hit')
    except:
        os.makedirs(pwd+'/cd-hit')
    os.chdir(pwd+'/cd-hit')
    if transdecoder == 'True' and len(transcripts) > 0:
        filtered_fastas = glob.glob(pwd+'/miniprot/*.filtered.fasta')+glob.glob(pwd+'/transdecoder/*.filtered.fasta')
    else:
        filtered_fastas = glob.glob(pwd+'/miniprot/*.filtered.fasta')
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
        spend = round(end_time - start_time)
        out_log.write('       '+filtered_fasta.rsplit('/')[-1]+': '+str(timedelta(seconds=spend))+'\n')
        out_log.flush()


out_log.write('\n#Orthofinder processing\n')
now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
start_time = time.time()
out_log.write(str(now)+'\n')
out_log.flush()
full_abbr = {}
try:
    size = os.path.getsize(pwd+'/'+prefix+'.'+str(num)+'.orthofinder')
    if len(glob.glob(pwd+'/'+prefix+'.'+str(num)+'.orthofinder/OrthoFinder/Results_*/Orthogroup_Sequences/*.fa')) < 10 or same == 'False':
        #os.system('rm -r '+pwd+'/'+prefix+'.'+str(num)+'.orthofinder')
        print('Please check the orthofinder results. We detect a older version. If you want to run a new one, just delete the '+pwd+'/'+prefix+'.'+str(num)+'.orthofinder')
        sys.exit()
        raise Exception
    else:
        with open(pwd+'/'+prefix+'.'+str(num)+'.orthofinder/Fullname_abbr.txt') as tmp:
            for line in tmp:
                tmp_n = line.strip('\n').split('\t')
                full_abbr[tmp_n[0]] = tmp_n[1]
    if len(full_abbr) != num:
        print('Please check the Full_name.txt in '+pwd+'/'+prefix+'.'+str(num)+'.orthofinder')
        sys.exit()
except:
    os.makedirs(pwd+'/'+prefix+'.'+str(num)+'.orthofinder')
    os.chdir(pwd+'/'+prefix+'.'+str(num)+'.orthofinder')
    cd_hit_fastas = glob.glob(pwd+'/cd-hit/*.cd-hit-0.85.fasta')+proteins
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
        cd_hit_fastas_1 = []
        for cd_hit_fasta in cd_hit_fastas:
            if cd_hit_fasta.rsplit('/')[-1].split('.')[0] in set(species_list):
                if 'transcript' not in cd_hit_fasta.rsplit('/')[-1]:
                    cd_hit_fastas_1.append(cd_hit_fasta)
                else:
                    if transdecoder == 'True' and 'transdecoder' in cd_hit_fasta.rsplit('/')[-1]:
                        cd_hit_fastas_1.append(cd_hit_fasta)
                    elif transdecoder == 'False' and 'transdecoder' not in cd_hit_fasta.rsplit('/')[-1]:
                        cd_hit_fastas_1.append(cd_hit_fasta)
        for cd_hit_fasta in cd_hit_fastas_1:
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
                            if '.' not in sequence and '*' not in sequence[:-1] and len(sequence) > seq_thres:
                                gene_name_1 = str.upper(sp_name)+'|Contig'+str(count).zfill(5)
                                out_1.write(gene_name+'\t'+gene_name_1+'\n')
                                out_2.write('>'+gene_name_1+'\n'+sequence+'\n')
                                count +=1
                out_1.close()
                out_2.close()
                num+=1
    os.chdir(pwd)
    if len(full_abbr.values()) != len(set(full_abbr.values())):
        print('It seems to have some duplicated ID that from the same_species. Please rename them and rerun.')
        sys.exit()
    orthofinder_cmd = 'orthofinder -f '+prefix+'.'+str(num)+'.orthofinder -og -t '+str(cores)
    os.system(orthofinder_cmd)
end_time = time.time()
spend = round(end_time - start_time)
if spend < 1:
    out_log.write('\n       existed and skipped\n')
else:
    out_log.write('\n       Running time: '+str(timedelta(seconds=spend))+'\n')
out_log.flush()


out_log.write('\n#Phylogenomic analysis\n')
now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
out_log.write(str(now)+'\n')
out_log.flush()
try:
    OGs_in = glob.glob(wd+'*.fa')
    if len(OGs_in ) < 10 or same == 'False':
        raise Exception
except:
    try:
        if os.path.getsize(wd) > 10:
            os.system('rm -r '+wd)
    except:
        print('We are trying to processing the trimming process for tree construction.')
    os.makedirs(wd)
    os.chdir(wd)
    OGs_wd = glob.glob(pwd+'/'+prefix+'.'+str(num)+'.orthofinder/OrthoFinder/Results_*/Orthogroup_Sequences/')[0]
    if len(glob.glob(pwd+'/'+prefix+'.'+str(num)+'.orthofinder/OrthoFinder/Results_*/Orthogroup_Sequences/*.fa')) < 10:
        ERROR = 'Orthofinder seems to not be finished properly, please check it'
        print(ERROR)
        out_log.write(ERROR)
        out_log.flush()
        sys.exit()
    OG_num = 0
    while OG_num < 40000:
        tmp_num = str(OG_num).zfill(7)
        OG_num += 1
        OG_name = OGs_wd+'/OG'+tmp_num+'.fa'
        with open(OG_name) as tmp_file:
            num_OGs = len(tmp_file.read().split('\n>'))
            if num_OGs < float(num*occupancy):
                break
            os.system('ln -s '+OG_name+' ./')


os.chdir(wd)
if check_point('check_occupancy_1st@@'+same+'@@'+str(round(occupancy,2))):
    out_log.write('       Check occupancy ('+str(round(occupancy*num))+'/'+str(num)+') 1st: existed and skipped.\n')
else:
    if os.path.isdir(wd+'rejected_few_taxa_1'):
        os.system('rm -r '+wd+'rejected_few_taxa_1')
    os.makedirs(wd+'rejected_few_taxa_1')
    if os.path.isdir(wd+'01.backup_all_OGs'):
        os.system('rm -r '+wd+'01.backup_all_OGs')
    os.makedirs(wd+'01.backup_all_OGs')
    start_time = time.time()
    OGs = glob.glob(wd+'*.fa')
    if len(OGs) < 10:
        print('Not enough files found in '+wd+'. Please check it')
    for OG in OGs:
        os.system(' cp '+OG+' '+wd+'01.backup_all_OGs')
        taxa = {}
        length = []
        with open(OG) as file:
            n = file.read().split('\n>')
            if len(n) >= float(num*occupancy):
                for i in n:
                    i = i.lstrip('>')
                    taxa_sp = i.split('|')[0]
                    seq = i.split('\n',1)[1].replace('\n','')
                    if taxa_sp not in taxa.keys():
                        taxa[taxa_sp] = 1
                    else:
                        taxa[taxa_sp] = taxa[taxa_sp] + 1
                    length.append(len(seq))
        if len(taxa) < float(num*occupancy)-0.5:
            os.system('mv '+OG+' '+wd+'rejected_few_taxa_1')
        else:
            remove_linebreak(OG)
    end_time = time.time()
    spend = round(end_time - start_time)
    if len(glob.glob('*.fa')) < 10:
        ERROR = '       Check occupancy ('+str(round(occupancy*num))+'/'+str(num)+') 1st: error! Please check it'
        print(ERROR)
        out_log.write(ERROR)
        out_log.flush()
        for backup in glob.glob(wd+'01.backup_all_OGs/*.fa'):
            os.system('mv '+backup+' ./')
        sys.exit()
    out_log.write('       Check occupancy ('+str(round(occupancy*num))+'/'+str(num)+') 1st: '+str(timedelta(seconds=spend))+'\n')
    write_check('check_occupancy_1st')
    out_log.flush()


if check_point('uniqHaplo@@'+same+'@@'+str(round(occupancy,2))):
    out_log.write('       Remove redundant sequences using uniqHaplo: existed and skipped.\n')
else:
    if os.path.isdir(wd+'02.backup_preUniqHaplo'):
        os.system('rm -r '+wd+'02.backup_preUniqHaplo')
    os.makedirs(wd+'02.backup_preUniqHaplo')
    start_time = time.time()
    def UniqHaplq(fa):
        os.system('cp '+fa+' ./02.backup_preUniqHaplo')
        os.system("perl "+uniqHaplo_path+" -a "+fa+' > '+fa+'.uniq')
        os.system('rm '+fa)
        os.system('mv '+fa+'.uniq '+fa)
    start_time = time.time()
    fas_preUniqHaplo = glob.glob(wd+'*.fa')
    with multiprocessing.Pool(processes=cores) as pool:
        pool.map(UniqHaplq, fas_preUniqHaplo)
    end_time = time.time()
    spend = round(end_time - start_time)
    if len(glob.glob('*.fa')) < 10:
        ERROR = '       Remove redundant sequences using uniqHaplo: error! Please check it'
        print(ERROR)
        out_log.write(ERROR)
        out_log.flush()
        for backup in glob.glob(wd+'02.backup_preUniqHaplo/*.fa'):
            os.system('mv '+backup+' ./')
        sys.exist()
    out_log.write('       Remove redundant sequences using uniqHaplo: '+str(timedelta(seconds=spend))+'\n')
    write_check('uniqHaplo')
    out_log.flush()


def mafft_aln(fa):
    os.system('cp '+fa+' ./03.backup_alignments')
    os.system("mafft --auto --localpair --quiet --maxiterate 1000 "+fa+' > '+fa+'.aln')
    os.system('mv '+fa+'.aln '+fa)
if check_point('Mafft@@'+same+'@@'+str(round(occupancy,2))):
    out_log.write('       Mafft alignment: existed and skipped.\n')
else:
    fas_mafft = glob.glob(wd+'*.fa')
    if os.path.isdir(wd+'03.backup_alignments'):
        os.system('rm -r '+wd+'03.backup_alignments')
    os.makedirs(wd+'03.backup_alignments')
    start_time = time.time()
    with multiprocessing.Pool(processes=cores) as pool:
        pool.map(mafft_aln, fas_mafft)
    end_time = time.time()
    spend = round(end_time - start_time)
    if len(glob.glob('*.fa')) < 10:
        ERROR = '       Mafft alignment: error! Please check it'
        print(ERROR)
        out_log.write(ERROR)
        out_log.flush()
        for backup in glob.glob(wd+'03.backup_alignments/*.fa'):
            os.system('mv '+backup+' ./')
        sys.exist()
    out_log.write('       Mafft alignment: '+str(timedelta(seconds=spend))+'\n')
    write_check('Mafft')
    out_log.flush()


def Hmmcleaner(fa):
    os.system('cp '+fa+' ./04.back_pre_HmmCleaner')
    os.system('HmmCleaner.pl '+fa+' --specificity')
    os.system('mv '+fa+' ./HmmCleaner_files')
    os.system('mv '+fa.replace('.fa','_hmm.log')+' ./HmmCleaner_files')
    os.system('mv '+fa.replace('.fa','_hmm.score')+' ./HmmCleaner_files')
    os.system('mv '+fa.replace('.fa','_hmm.fasta')+' '+fa)
if check_program_exist('HmmCleaner.pl'):
    if check_point('HmmCleaner.pl@@'+same+'@@'+str(round(occupancy,2))): 
        out_log.write('       Clean alignments with HmmCleaner: existed and skipped.\n')
    else:
        fas_hmmcleaner = glob.glob(wd+'*.fa')
        if os.path.isdir(wd+'04.back_pre_HmmCleaner'):
            os.system('rm -r '+wd+'04.back_pre_HmmCleaner')
        os.makedirs(wd+'04.back_pre_HmmCleaner')
        if os.path.isdir(wd+'HmmCleaner_files'):
            os.system('rm -r '+wd+'HmmCleaner_files')
        os.makedirs(wd+'HmmCleaner_files')
        start_time = time.time()
        with multiprocessing.Pool(processes=cores) as pool:
            pool.map(Hmmcleaner, fas_hmmcleaner)
        end_time = time.time()
        spend = round(end_time - start_time)
        if len(glob.glob('*.fa')) < 10:
            ERROR = '       Clean alignments with HmmCleaner: error! Please check it'
            print(ERROR)
            out_log.write(ERROR)
            out_log.flush()
            for backup in glob.glob(wd+'04.back_pre_HmmCleaner/*.fa'):
                os.system('mv '+backup+' ./')
            sys.exist()
        out_log.write('       Clean alignments with HmmCleaner: '+str(timedelta(seconds=spend))+'\n')
        write_check('HmmCleaner.pl')
else:
    out_log.write('       Clean alignments with HmmCleaner: HmmCleaner Not Found and skipped.\n')
out_log.flush()    


def trimal(fa):
    os.system('cp '+fa+' ./05.backup_pre-trimal')
    os.system('trimal -in '+fa+' -out '+fa+'.trimal -automated1')
    os.system('mv '+fa+'.trimal '+fa)
if check_point('trimal@@'+same+'@@'+str(round(occupancy,2))):
    out_log.write('       Trim alignments with trimal: existed and skipped.\n')
else:
    fas_trimal = glob.glob(wd+'*.fa')
    if os.path.isdir(wd+wd+'05.backup_pre-trimal'):
        os.system('rm -rf '+wd+'05.backup_pre-trimal')
    os.makedirs(wd+'05.backup_pre-trimal')
    start_time = time.time()
    with multiprocessing.Pool(processes=cores) as pool:
        pool.map(trimal, fas_trimal)
    end_time = time.time()
    spend = round(end_time - start_time)
    if len(glob.glob('*.fa')) < 10:
        ERROR = '       Trim alignments with trimal: error! Please check it'
        print(ERROR)
        out_log.write(ERROR)
        out_log.flush()
        for backup in glob.glob(wd+'05.backup_pre-trimal/*.fa'):
            os.system('mv '+backup+' ./')
        sys.exist()
    out_log.write('       Trim alignments with trimal: '+str(timedelta(seconds=spend))+'\n')
    write_check('trimal')
    out_log.flush()



def BMGE(fa):
    os.system('cp '+fa+' ./06.backup_pre-BMGE')
    #os.system('java -jar '+BMGE_path+' -t AA -i '+fa+' -of '+fa+'.BMGE')
    os.system('bmge -t AA -i '+fa+' -of '+fa+'.BMGE')
    os.system('mv '+fa+'.BMGE '+fa)
if check_point('BMGE@@'+same+'@@'+str(round(occupancy,2))):
    out_log.write('       Trim alignments with BMGE: existed and skipped.\n')
else:
    fas_BMGE = glob.glob(wd+'*.fa')
    if os.path.isdir(wd+wd+'06.backup_pre-BMGE'):
        os.system('rm -rf '+wd+'06.backup_pre-BMGE')
    os.makedirs(wd+'06.backup_pre-BMGE')
    start_time = time.time()
    with multiprocessing.Pool(processes=cores) as pool:
        pool.map(BMGE, fas_BMGE)
    end_time = time.time()
    spend = round(end_time - start_time)
    if len(glob.glob('*.fa')) < 10:
        ERROR = '       Trim alignments with BMGE: error! Please check it'
        print(ERROR)
        out_log.write(ERROR)
        out_log.flush()
        for backup in glob.glob(wd+'06.backup_pre-BMGE/*.fa'):
            os.system('mv '+backup+' ./')
        sys.exist()
    out_log.write('       Trim alignments with BMGE: '+str(timedelta(seconds=spend))+'\n')
    write_check('BMGE')
    out_log.flush()

def AlignmentCompare(fa):
    os.system('cp '+fa+' ./07.back_pre_AlignmentCompare')
    os.system('java -cp '+AlignmentCompare_path+' AlignmentCompare '+fa)
if check_point('AlignmentCompare@@'+same+'@@'+str(round(occupancy,2))):
    out_log.write('       AlignmentCompare (<20 amino acids): existed and skipped.\n')
else:
    fas_aligncompare = glob.glob(wd+'*.fa')
    if os.path.isdir(wd+'07.back_pre_AlignmentCompare'):
        os.system('rm -r '+wd+'07.back_pre_AlignmentCompare')
    os.makedirs(wd+'07.back_pre_AlignmentCompare')
    start_time = time.time()
    with multiprocessing.Pool(processes=cores) as pool:
        pool.map(AlignmentCompare, fas_aligncompare)
    end_time = time.time()
    spend = round(end_time - start_time)
    if len(glob.glob('*.fa')) < 10:
        ERROR = '       AlignmentCompare (<20 amino acids): error! Please check it'
        print(ERROR)
        out_log.write(ERROR)
        out_log.flush()
        for backup in glob.glob(wd+'07.back_pre_AlignmentCompare/*.fa'):
            os.system('mv '+backup+' ./')
        sys.exit()
    out_log.write('       AlignmentCompare (<20 amino acids): '+str(timedelta(seconds=spend))+'\n')
    write_check('AlignmentCompare')
    out_log.flush()

if check_point('check_occupancy_2nd@@'+same+'@@'+str(round(occupancy,2))):
    out_log.write('       Check occupancy ('+str(round(occupancy*num))+'/'+str(num)+') 2nd: existed and skipped.\n')
else:
    if os.path.isdir(wd+'rejected_few_taxa_2'):
        os.system('rm -r '+wd+'rejected_few_taxa_2')
    os.makedirs(wd+'rejected_few_taxa_2')
    if os.path.isdir(wd+'08.backup_check_occupancy_2nd'):
        os.system('rm -r '+wd+'08.backup_check_occupancy_2nd')
    os.makedirs(wd+'08.backup_check_occupancy_2nd')
    start_time = time.time()
    OGs = glob.glob(wd+'*.fa')
    for OG in OGs:
        os.system('cp '+OG+' ./08.backup_check_occupancy_2nd')
        taxa = {}
        length = []
        if os.path.getsize(OG) < 100:
            os.system('mv '+OG+' '+wd+'rejected_few_taxa_2')
        else:
            try:
                with open(OG) as file:
                    n = file.read().split('\n>')
                    if len(n) >= float(num*occupancy):
                        for i in n:
                            i = i.lstrip('>')
                            taxa_sp = i.split('|')[0]
                            seq = i.split('\n',1)[1].replace('\n','')
                            if taxa_sp not in taxa.keys():
                                taxa[taxa_sp] = 1
                            else:
                                taxa[taxa_sp] = taxa[taxa_sp] + 1
                            length.append(len(seq))
                if len(taxa) < float(num*occupancy)-0.5 or len(seq) < seq_thres:
                    os.system('mv '+OG+' '+wd+'rejected_few_taxa_2')
                else:
                    remove_linebreak(OG)
            except:
                os.system('mv '+OG+' '+wd+'rejected_few_taxa_2')
                continue
    end_time = time.time()
    spend = round(end_time - start_time)
    if len(glob.glob('*.fa')) < 10:
        ERROR = '       Check occupancy ('+str(round(occupancy*num))+'/'+str(num)+') 2nd: error! Please check it'
        print(ERROR)
        out_log.write(ERROR)
        out_log.flush()
        for backup in glob.glob(wd+'07.backup_check_occupancy_2nd/*.fa'):
            os.system('mv '+backup+' ./')
        sys.exit()
    out_log.write('       Check occupancy ('+str(round(occupancy*num))+'/'+str(num)+') 2nd: '+str(timedelta(seconds=spend))+'\n')
    write_check('check_occupancy_2nd')
    out_log.flush()

def iqtree2(fa):
    iqtree2_cmd = 'iqtree2 -B 1000 -T 1 -s '+fa
    tre = fa.replace('.fa','.tre')
    contree = fa+'.contree'
    if os.path.isfile(tre) is False:
        os.system(iqtree2_cmd)
        os.system('mv '+contee+' '+tre)
    elif os.path.isfile(contree) is False:
        os.system('mv '+contee+' '+tre)
        

start_time = time.time()
fas_tree = glob.glob(wd+'*.fa')
if MODE == 'FastTree':
    #FastTree for each OG
    fas_tree = glob.glob(wd+'*.fa')
    out_fasttree = open('run_fasttree.sh','w')
    out_fasttree.write('export OMP_NUM_THREADS='+str(cores)+'\n')
    for fa_fasttree in fas_tree:
        tre_name = fa_fasttree.rsplit('.fa',1)[0]+'.tre'
        if os.path.isfile(tre_name) is False:
            out_fasttree.write('FastTreeMP -slow -gamma '+ fa_fasttree +' > '+fa_fasttree+'.tre\n')
            out_fasttree.write('mv '+fa_fasttree+'.tre '+tre_name+'\n')
        else:
            if os.path.getsize(tre_name) > 10:
                out_fasttree.write('#FastTreeMP -slow -gamma '+ fa_fasttree +' > '+fa_fasttree+'.tre\n')
                out_fasttree.write('#mv '+fa_fasttree+'.tre '+tre_name+'\n')
    out_fasttree.close()
    os.system('sh run_fasttree.sh')
else:
    with multiprocessing.Pool(processes=cores) as pool:
        pool.map(iqtree2, fas_tree)

end_time = time.time()
spend = round(end_time - start_time)
if spend < 1:
    out_log.write('       Gene Tree in each OG ('+MODE+'): existed and skipped\n')
else:
    out_log.write('       Gene Tree in each OG ('+MODE+'):'+str(timedelta(seconds=spend))+'\n')
out_log.flush()

fas = glob.glob(wd+'*.fa')
for fa in fas:
    length = []
    with open(fa) as tmp_file:
        tmp_n = tmp_file.read().split('\n>')
        if len(tmp_n) > round(occupancy,2)*num-0.5:
            for tmp_i in tmp_n:
                if '|' not in tmp_i:
                    os.system('rm '+fa)
                    os.system('rm '+fa.rsplit('.',1)[0]+'.tre')
                else:
                    length.append(len(tmp_i.split('\n',1)[1].replace('\n','')))
        if len(set(length)) !=1 or len(tmp_i.split('\n',1)[1].replace('\n','')) < seq_thres:
            os.system('rm '+fa)
            os.system('rm '+fa.rsplit('.',1)[0]+'.tre')

start_time = time.time()
phylopypruner_cmd = 'phylopypruner --threads '+str(cores)+' --min-taxa '+ str(round(occupancy*num))+' --min-len '+str(seq_thres)+' --dir . --min-support 0.75 --mask pdist --trim-divergent 0.75 --min-pdist 0.01 --prune MI'
#out_log.write('       phylopypruner: '+phylopypruner_cmd+'\n')
out_log.flush()
try:
    tmp_test = glob.glob(wd+'phylopypruner_output/filtered/*.fa')
    if len(tmp_test) < 10 or same == 'False':
        raise Exception
except:
    os.system(phylopypruner_cmd)
    end_time = time.time()
    spend = round(end_time - start_time)
    out_log.write('       phylopypruner: '+str(timedelta(seconds=spend))+'\n')
    out_log.flush()

    #check duplicated OG and output non-redudant OGs.
    matrix = {}
    new_matrix = {}
    os.makedirs(wd+'phylopypruner_output/filtered')
    out_fas = open(wd+'phylopypruner_output/supermatrix.new.fas','w')
    out_partition = open(wd+'phylopypruner_output/partition_data.new.txt','w')
    with open(wd+'phylopypruner_output/supermatrix.fas') as file:
        content = file.read().lstrip('>').split('\n>')
        for i in content:
            matrix_name = i.split('\n',1)[0]
            matrix_seq = i.split('\n',1)[1].replace('\n','')
            matrix[matrix_name] = matrix_seq

    partitions = {}
    duplicated = 0
    with open(wd+'phylopypruner_output/partition_data.txt') as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip('\n')
            OG_name = line.split(' ')[1]
            locations = line.split(' ')[-1].split('-')
            new_start = str(int(locations[0]) - duplicated)
            new_end = str(int(locations[1]) - duplicated)
            new_locations = new_start+'-'+new_end
            if OG_name not in partitions.keys():
                partitions[OG_name] = new_locations
                out_partition.write('AUTO, '+OG_name+' = '+new_locations+'\n')
                tmp_OG = open(wd+'phylopypruner_output/filtered/'+OG_name+'_pruned.fa','w')
                for species in matrix.keys():
                    if species not in new_matrix.keys():
                        new_matrix[species] = matrix[species][int(locations[0])-1:int(locations[1])]
                    else:
                        new_matrix[species] = new_matrix[species] + matrix[species][int(locations[0])-1:int(locations[1])]
                    count_X = matrix[species][int(locations[0])-1:int(locations[1])].count('X')+matrix[species][int(locations[0])-1:int(locations[1])].count('-')+matrix[species][int(locations[0])-1:int(locations[1])].count('?')
                    if count_X < len(matrix[species][int(locations[0])-1:int(locations[1])]):
                        tmp_OG.write('>'+full_abbr[species]+'\n'+matrix[species][int(locations[0])-1:int(locations[1])]+'\n')
                tmp_OG.close()
            else:
                duplicated = duplicated + 1 + int(locations[1]) - int(locations[0])
        print(duplicated)
        for species in new_matrix.keys():
            out_fas.write('>'+full_abbr[species]+'\n'+new_matrix[species]+'\n')
    out_fas.close()
    out_partition.close()

def fasta2phy(f):
    longest_name = ''
    taxa = 0
    site = 0
    with open(f) as file:
        out = open(f.rsplit('.',1)[0]+'.phy','w')
        n = file.read().lstrip('>').strip('\n').split('\n>')
        for i in n:
            name = i.split('\n',1)[0]
            seq = i.split('\n',1)[1].replace('\n','')
            if len(name) >= len(longest_name):
                longest_name = name
            taxa += 1
        site = len(seq)
        gap = longest_name+'aaaa'
        out.write(str(taxa)+' '+str(site)+'\n')
        for i in n:
            name = i.split('\n',1)[0]
            seq = i.split('\n',1)[1].replace('\n','')
            out.write(name.ljust(len(gap))+seq+'\n')
        out.close()

OGs = glob.glob(wd+'phylopypruner_output/filtered/*_pruned.fa')
random.shuffle(OGs)
sample_size = 0
sub_1 = {}
out_sub1_fa = open(wd+'phylopypruner_output/'+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.2500000.fa','w')
out_sub1_txt = open(wd+'phylopypruner_output/'+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.2500000.txt','w')
out_sub2_fa = open(wd+'phylopypruner_output/'+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.5000000.fa','w')
out_sub2_txt = open(wd+'phylopypruner_output/'+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.5000000.txt','w')
out_sub1_start = 1
out_sub2_start = 1
sub_2 = {}
for species in full_abbr.values():
    longest_sp = ''
    if len(species) >= len(longest_sp):
        longest_sp = species
for OG in OGs:
    with open(OG) as tmp_file:
        tmp_dict = {}
        OG_name = OG.split('/')[-1].split('_pruned.fa')[0]
        tmp_n = tmp_file.read().lstrip('>').strip('\n').split('\n>')
        for tmp_i in tmp_n:
            tmp_seq = tmp_i.split('\n',1)[1].replace('\n','')
            tmp_name = tmp_i.split('\n',1)[0]
            tmp_dict[tmp_name] = tmp_seq
            tmp_length = len(tmp_seq)
        if sample_size < 2500000:
            out_sub1_end = out_sub1_start+tmp_length-1
            new_locations1 = str(out_sub1_start)+' - '+str(out_sub1_end)
            out_sub1_txt.write('AUTO, '+OG_name+' = '+new_locations1+'\n')
            for sp in full_abbr.values():
                if sp not in sub_1.keys():
                    if sp in tmp_dict:
                        sub_1[sp] = tmp_dict[sp]
                    else:
                        sub_1[sp] = 'X'*tmp_length
                else:
                    if sp in tmp_dict:
                        sub_1[sp] = sub_1[sp]+tmp_dict[sp]
                    else:
                        sub_1[sp] = sub_1[sp]+'X'*tmp_length
            out_sub1_start += tmp_length
        if sample_size < 5000000:
            out_sub2_end = out_sub2_start+tmp_length-1
            new_locations2 = str(out_sub2_start)+' - '+str(out_sub2_end)
            out_sub2_txt.write('AUTO, '+OG_name+' = '+new_locations2+'\n')
            for sp in full_abbr.values():
                if sp not in sub_2.keys():
                    if sp in tmp_dict.keys():
                        sub_2[sp] = tmp_dict[sp]
                    else:
                        sub_2[sp] = 'X'*tmp_length
                else:
                    if sp in tmp_dict.keys():
                        sub_2[sp] = sub_2[sp]+tmp_dict[sp]
                    else:
                        sub_2[sp] = sub_2[sp]+'X'*tmp_length
            out_sub2_start += tmp_length
        sample_size += num*tmp_length
for species in full_abbr.values():
    out_sub1_fa.write('>'+species+'\n'+sub_1[species]+'\n')
    out_sub2_fa.write('>'+species+'\n'+sub_2[species]+'\n')
out_sub1_txt.close()
out_sub2_txt.close()
out_sub1_fa.close()
out_sub2_fa.close()
fasta2phy(wd+'phylopypruner_output/'+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.2500000.fa')
fasta2phy(wd+'phylopypruner_output/'+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.5000000.fa')

with open(wd+'phylopypruner_output/run_phylobayes.2500000.sh','w') as tmp_out:
    tmp_out.write('#Three chains running in Phylobayes\n')
    for tmp in range(1,4):
        #mpirun -np '+str(thread)+' pb_mpi -dc -cat -gtr -dgam 4 -d '+input_file+' '+prefix+'.PB.chain'+str(i)
        tmp_out.write('mpirun -np 10 pb_mpi -dc -cat -gtr -dgam 4 -d '+wd+'phylopypruner_output/'+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.2500000.phy '+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.2500000.PB.chain'+str(tmp)+'\n')
    tmp_out.write('#Check maxdiff: if maxdiff < 0.1, it would be a good run\n')
    tmp_out.write('#bpcomp -x 300 '+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.2500000.PB.chain*.treelist\n')
with open(wd+'phylopypruner_output/run_phylobayes.5000000.sh','w') as tmp_out:
    tmp_out.write('#Three chains running in Phylobayes\n')
    for tmp in range(1,4):
        tmp_out.write('mpirun -np 10 pb_mpi -dc -cat -gtr -dgam 4 -d '+wd+'phylopypruner_output/'+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.2500000.phy '+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.5000000.PB.chain'+str(tmp)+'\n')
    tmp_out.write('#Check maxdiff: if maxdiff < 0.1, it would be a good run\n')
    tmp_out.write('#bpcomp -x 300 '+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.2500000.PB.chain*.treelist\n')

start_time = time.time()
os.chdir(wd+'phylopypruner_output')
fasttree_cmd = 'FastTreeMP -slow -gamma supermatrix.new.fas > FastTree.tre'
out_log.flush()
try:
    if os.path.getsize(wd+'phylopypruner_output/'+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.FastTree.full.tre') < 10:
        raise Exception
except:
    with open('run_fasttree.sh','w') as tmp_out:
        tmp_out.write('export OMP_NUM_THREADS='+str(cores)+'\n')
        tmp_out.write(fasttree_cmd)
    os.system('sh run_fasttree.sh')

os.system('cp '+wd+'phylopypruner_output/FastTree.tre '+wd+'phylopypruner_output/'+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.FastTree.full.tre')
os.system('ln -fs '+wd+'phylopypruner_output/'+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.FastTree.full.tre '+pwd+'/'+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.FastTree.full.tre')
end_time = time.time()
spend = round(end_time - start_time)
if spend < 1:
    out_log.write('       FastTree for supermatrix.fas: existed and skipped\n')
else:
    out_log.write('       FastTree for supermatrix.fas: '+str(timedelta(seconds=spend))+'\n')
out_log.flush()


start_time = time.time()
iqtree2_cmd = 'iqtree2 --threads-max '+str(cores)+' -T AUTO -B 1000 -s supermatrix.new.fas -Q partition_data.new.txt -m MFP'
#out_log.write('       IQ-TREE2 -m MFP construction: '+iqtree2_cmd+'\n')
out_log.flush()
try:
    if os.path.getsize(wd+'phylopypruner_output/'+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.IQTREE2.full.tre') < 10:
        raise Exception
except:
    os.system(iqtree2_cmd)

os.system('cp '+wd+'phylopypruner_output/partition_data.new.txt.contree '+wd+'phylopypruner_output/'+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.IQTREE2.full.tre')
os.system('ln -fs '+wd+'phylopypruner_output/'+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.IQTREE2.full.tre '+pwd+'/'+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.IQTREE2.full.tre')
end_time = time.time()
spend = round(end_time - start_time)
if spend < 1:
    out_log.write('       IQ-TREE2 -m MFP construction: existed and skipped.\n')
else:
    out_log.write('       IQ-TREE2 -m MFP construction: '+str(timedelta(seconds=spend))+'\n')
out_log.flush()

now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
out_log.write('\nAll done: '+str(now)+'\n')
et = time.time()
spend = round(et - st,1)
out_log.write('\nSpend time: '+str(timedelta(seconds=spend))+'\n')


with open('run_ASTRAL_FastTreeMP.sh','w') as tmp_out:
    tmp_out.write('export OMP_NUM_THREADS='+str(cores)+'\n')
    for OG in OGs:
        tre_name = OG+'.tre'
        if os.path.isfile(tre_name) is False:
            tmp_out.write('FastTreeMP -slow -gamma '+ OG +' > '+tre_name+'\n')
        elif os.path.getsize(tre_name) < 10:
            tmp_out.write('FastTreeMP -slow -gamma '+ OG +' > '+tre_name+'\n')
        else:
            tmp_out.write('#FastTreeMP -slow -gamma '+ OG +' > '+tre_name+'\n')
    tmp_out.write('cat '+wd+'phylopypruner_output/filtered/*.tre > '+wd+'phylopypruner_output/'+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.genetree_FastTreeMP.tre\n')
    tmp_out.write('astral -i '+wd+'phylopypruner_output/'+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.genetree_FastTreeMP.tre'+' -o '+wd+'phylopypruner_output/'+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+'.astral_FastTreeMP\n')


with open('run_ASTRAL_IQTREE2.sh','w') as tmp_out:
    tmp_out.write('python3 '+path+'/scripts/run_ASTRAL_IQTREE2.py '+prefix+'.'+str(num)+'__'+str(round(occupancy,2))+' '+str(cores))

print('\nAll done: '+str(now)+'\n')
print('Spend time: '+str(timedelta(seconds=spend))+'\n')
out_log.write('\nPlease check the constructed trees: '+prefix+'.FastTree.full.tre and '+prefix+'.IQTREE2.full.tre')
print('Please check the constructed trees: '+prefix+'.FastTree.full.tre and '+prefix+'.IQTREE2.full.tre')
print('And commands for ASTRAL (two ways: FastTreeMP and IQTREE2) and Phylobayes tree inference could be accessed in '+wd+'/phylopypruner_output')
out_log.close()

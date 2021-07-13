#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import shutil
import re
pd.set_option('display.max_colwidth', None)


#func with filtering contigs by len and cov
def check_graph(id_name,table,new_dir):
    two = table[table['ID'] == id_name]
    two = two.reset_index()
    two['cov_log'] = np.log10(two['COV'])
    two = two.sort_values('cov_log',ascending = True)
    two['diff']=two['cov_log'] -two['cov_log'].shift(1)
    two = two.fillna(0)
    two = two.sort_values('diff', ascending = True)
    limit = two.loc[two['diff'].idxmax(),'COV']
    limit_log = two.loc[two['diff'].idxmax(),'cov_log']
    plt.scatter(x = two['LENGTH'], y = two['cov_log'])
    id_true_name = int(two.iloc[0]['ID'])
    plt.title(str(id_true_name) + ' ' + 'limit=' + str(round(limit, 2)) + ' ' + 'log_limit=' + str(round(limit_log, 2)))
    plt.savefig( graphs_before_folder +'/' +'%s.png'% str(two.iloc[0]['ID']))
    #plt.show()
    if limit < 100:
        limit_data[id_name] = limit
    return two

#if sum len is ok for this type
def check_length(data_new, len_col,type_col,len_alarm_col):
    data_new[len_alarm_col] = 0
    data_new.loc[(data_new[type_col] == 'Pseudomonas aeruginosa') & ((data_new[len_col] < 5500000) |(data_new[len_col] > 7500000)), len_alarm_col] = 1
    data_new.loc[(data_new[type_col] == 'Acinetobacter baumannii') & ((data_new[len_col] < 3500000) |(data_new[len_col] > 5000000)), len_alarm_col] = 1
    data_new.loc[(data_new[type_col] == 'Klebsiella pneumoniae') & ((data_new[len_col] < 5100000) |(data_new[len_col] > 6600000)), len_alarm_col] = 1
    data_new.loc[(data_new[type_col] == 'Escherichia coli') & ((data_new[len_col] < 4600000) |(data_new[len_col] > 6100000)), len_alarm_col] = 1
    data_new.loc[(data_new[type_col] == 'Serratia marcescens') & ((data_new[len_col] < 4700000) |(data_new[len_col] > 6200000)), len_alarm_col] = 1
    data_new.loc[(data_new[type_col] == 'Enterobacter cloacae') & ((data_new[len_col] < 4500000) |(data_new[len_col] > 6000000)), len_alarm_col] = 1
    data_new.loc[(data_new[type_col] == 'Staphylococcus aureus') & ((data_new[len_col] < 2300000) |(data_new[len_col] > 3800000)), len_alarm_col] = 1
    data_new.loc[(data_new[type_col] == 'Enterococcus faecium') & ((data_new[len_col] < 2500000) |(data_new[len_col] > 3900000)), len_alarm_col] = 1
    return data_new

#parse id column
def change_name(table, column):
    for index, row in table.iterrows():
        new_name = int(re.search(r'\d+', row[column]).group(0))
        table.at[index, column] = new_name
    return table

#cov and len information from contigs
def stat_from_fasta (path_genome, name_node, name_len, name_cov):
    contig = {}
    for path,dirs,files in os.walk(path_genome):
        for file in files:
            if file.startswith('.'):
                continue
            id_name = int(re.search(r'\d+',file).group(0))
            contig[id_name] = {}
            f =  open(path_genome + '/' + file, 'r',encoding="utf8", errors='ignore')
            for line in f:
                if line.startswith('>') == True:
                    line = line.rstrip('\n')
                    new = line.split('_')
                    contig[id_name][int(new[1])] = {}
                    contig[id_name][int(new[1])]['len'] = int(new[3])
                    contig[id_name][int(new[1])]['cov'] = float(new[5])
    contig_info = pd.DataFrame([(k,k1,v1['len'], v1['cov']) for k,v in contig.items() for k1,v1 in v.items()], columns = ['file','NODE','LENGTH', 'COV'])
    contig_info.columns = ['ID', name_node, name_len, name_cov]
    table = contig_info.groupby('ID').agg({name_len :'sum', name_cov : 'mean', name_node : 'count'}).reset_index()
    return contig_info,table

if( len( sys.argv ) != 4 ):
    print('Parse tables after TORMES pipeline and get summary table')
    print('Usage: python3 parse_tormes_results.py <path_folder_with_all_tormes_folders> <path_table_with_metadata> <output_path_sum_table.csv>')
    exit()

#download all tables with results. folder with all tormes results and table with metadata are from cmd
tormes_folder = sys.argv[1] #folder with tormes_results
genomes = tormes_folder + 'genomes'
graphs_before_folder = os.path.join(tormes_folder, 'cov_len_graphs_before')
os.mkdir(graphs_before_folder)
genomes_after_filter = os.path.join(tormes_folder, 'genomes_after_filter')
os.mkdir(genomes_after_filter)
meta = sys.argv[2]
metadata = pd.read_csv(meta,sep = '\t') # folder with metadata
metadata.columns = ['ID', 'type_from_smol', 'Город', 'Материал', 'Маркеры резистентности',
       'Concentration, ng/mkl', 'A260', '260/280', '260/230']
metadata_short = metadata[['ID', 'type_from_smol', 'Маркеры резистентности']]
mlst = pd.read_csv(tormes_folder + 'mlst/mlst.tab', sep = '\t', header = None)
mlst.columns = ['ID', 'type_mlst', 'st',3, 4, 5, 6, 7, 8, 9 ]
mlst_short = mlst[['ID', 'type_mlst', 'st']]
mlst_short = change_name(mlst_short, 'ID')
data_tax = pd.read_csv(tormes_folder + 'taxonomic_identification/taxonomic-identification-16S-rRNA.RDP.txt', sep = '\t')
data_tax = data_tax.groupby('Gene')['Genus'].apply(lambda x: "[%s]" % ', '.join(x)).reset_index()
data_tax = change_name(data_tax, 'Gene')
data_tax.columns = ['ID', 'Genus']

#do cat of all tab with amr_genes all samples
#only resfinder folder, should change to path /antibiotic_resistance_genes/resfinder/ - done
cmd = 'cat %s'%tormes_folder    + 'antibiotic_resistance_genes/resfinder/*.tab > %s'%tormes_folder+ 'antibiotic_resistance_genes/resfinder/all_genes.tab'
os.system(cmd)
amr = pd.read_csv(tormes_folder + 'antibiotic_resistance_genes/resfinder/all_genes.tab', sep = '\t')
amr = amr[amr['#FILE'] != '#FILE']
amr = change_name(amr, '#FILE')
amr.columns = ['ID', 'SEQUENCE', 'START', 'END', 'STRAND', 'GENE', 'COVERAGE',
       'COVERAGE_MAP', 'GAPS', '%COVERAGE', '%IDENTITY', 'DATABASE',
       'ACCESSION', 'PRODUCT', 'RESISTANCE']
amr_short = amr.groupby('ID')['PRODUCT'].apply(lambda x: ', '.join(x)).reset_index() #all amr_genes in list

contig_info,sum_table = stat_from_fasta(genomes, 'NODE', 'LENGTH', 'COV')

#limit graphs
sample_id = sum_table['ID'].tolist()
limit_data = {}
for name in sample_id:
    check_graph(name,contig_info,'cov_len_graphs_before' )
sum_table['limit'] = 0
for index, row in sum_table.iterrows():
    if row['ID'] in limit_data.keys():
        sum_table.loc[index, 'limit'] = limit_data[row['ID']]

#write new fasta with new limit of cov
for path,dirs,files in os.walk(genomes):
    for file in files:
        ext = os.path.splitext(file)[-1].lower()
        if ext == '.fasta':
            id_name = int(re.search(r'\d+',file).group(0))
            if id_name in limit_data.keys() :
                f =  open(genomes + '/' + file, 'r', encoding="utf8", errors = 'ignore')
                d = open(genomes_after_filter + '/' + "new_%s"%file[3:], "w")
                file_split = f.read().split('>')
                for node in file_split:
                    first = node.split('\n')
                    if len(first) > 1:
                        second = first[0].split('_')
                        if float(second[5]) >= limit_data[id_name]:
                            d.write('>' + node)
                d.close()
            else:
                source = genomes + '/' + file
                dest = genomes_after_filter + '/' + "new_%s"%file[3:]
                res = shutil.copy(source,dest)

#stat from new fasta
contig_info2,sum_table_new = stat_from_fasta(genomes_after_filter,'NEW_NODE', 'NEW_LENGTH', 'NEW_COV')

sum_table = pd.merge(sum_table, sum_table_new, on = 'ID')
sum_table = pd.merge(sum_table, metadata_short, on = 'ID')
sum_table = pd.merge(sum_table, mlst_short, on = 'ID')
sum_table = pd.merge(sum_table, data_tax, on = 'ID')

#flag if 16s is diff with meta
sum_table['flag_16s'] = 0
for index, row in sum_table.iterrows():
    genus = row['Genus'][1:-1].split(',')
    genus2 = row['type_from_smol'].split(' ')[0]
    flag = 0
    for name in genus:
        name = name.lstrip(' ')
        if name.startswith('unclassified'):
            continue
        elif name != genus2:
            flag+=1
    sum_table.at[index, 'flag_16s'] = flag

#flag if len is out of ranges
sum_table['len_alarm'] = 0
sum_table = check_length(sum_table, 'LENGTH', 'type_from_smol', 'len_alarm')

sum_table = pd.merge(sum_table, amr_short, on = 'ID')
sum_table.to_csv(sys.argv[3], index= False)

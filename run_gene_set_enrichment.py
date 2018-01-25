import numpy as np 
import os
import sys
import pdb
import gzip

def extract_ensamble_ids_from_qtl_results(qtl_results_file):
    f = open(qtl_results_file)
    head_count = 0
    hits = {}
    background = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        pvalue = float(data[-3])
        ensamble_id = data[5]
        if pvalue < .07:
            hits[ensamble_id] = 1
        else:
            background[ensamble_id] = 1
    return hits, background

def convert_from_ensamble_to_gene_symbol(ensamble_hits, gencode_file):
    f = gzip.open(gencode_file)
    gene_symbol_hits = []
    for line in f:
        line = line.rstrip()
        data = line.split()
        if line.startswith('#'):
            continue
        line_ensamble_id = data[9].split('"')[1]
        line_gene_symbol = data[17].split('"')[1]
        if line_ensamble_id in ensamble_hits:
            gene_symbol_hits.append(line_gene_symbol)
    return np.unique(gene_symbol_hits)

def sort_gsea(save_file, new_save_file):
    f = open(save_file)
    t = open(new_save_file,'w')
    pairs = []
    for i,line in enumerate(f):
        line = line.rstrip()
        data = line.split()
        if i < 4:
            continue
        pvalue = float(data[6])
        pairs.append((pvalue, line))
    sorted_pairs = sorted(pairs, key=lambda x: x[0])
    for pair in sorted_pairs:
        liner = pair[1]
        t.write(liner + '\n')
    t.close()

qtl_results_dir = sys.argv[1]
qtl_gene_set_enrichment_dir = sys.argv[2]
parameter_string = sys.argv[3]

gencode_file = '/project2/gilad/bstrober/ipsc_differentiation/preprocess_input_data/gencode.v19.annotation.gtf.gz'

qtl_results_file = qtl_results_dir + parameter_string + '_dynamic_qtl_results.txt'

ensamble_hits, ensamble_background = extract_ensamble_ids_from_qtl_results(qtl_results_file)

gene_symbol_hits = convert_from_ensamble_to_gene_symbol(ensamble_hits, gencode_file)

gene_symbol_background = convert_from_ensamble_to_gene_symbol(ensamble_background, gencode_file)


hits_file = qtl_gene_set_enrichment_dir + parameter_string + '_hit_genes.txt'

background_file = qtl_gene_set_enrichment_dir + parameter_string + '_background_genes.txt'


np.savetxt(hits_file, gene_symbol_hits,fmt="%s",delimiter="\n")
np.savetxt(background_file, gene_symbol_background,fmt="%s",delimiter="\n")


save_file = qtl_gene_set_enrichment_dir + parameter_string + '_gsea_output.txt'
geneset_file = '/project2/gilad/bstrober/tools/tools/gsea/data/' + 'c2.cp.kegg.v5.1.symbols.gmt.txt'
os.system('gsea ' + hits_file + ' ' + background_file + ' ' + geneset_file + ' ' + save_file)


new_save_file = qtl_gene_set_enrichment_dir + parameter_string + '_gsea_sorted_output.txt'
sort_gsea(save_file, new_save_file)
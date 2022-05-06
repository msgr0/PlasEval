#!/usr/bin/python
#
from __future__ import division
import os
import subprocess
import argparse

##USAGE: python analyse_sample.py path_to_PlasEval sample_id path_to_out_dir results_file
    
def evaluate(sample, out_dir, results, compare_script, extract_script, sample_script):
    p = subprocess.Popen('bash %s chr %s' % (sample_script, sample), stdout = subprocess.PIPE, shell = True)
    output, _ = p.communicate()
    p.wait()
    chromosomes = output.rstrip().decode().split(',')    

    p = subprocess.Popen('bash %s pla %s' % (sample_script, sample), stdout = subprocess.PIPE, shell = True)
    output, _ = p.communicate()
    p.wait()
    plasmids = output.rstrip().decode().split(',')    
    
    eval_dir = out_dir + '/eval'
    subprocess.call('mkdir %s' % eval_dir, shell = True)
    chr_references_fasta = eval_dir + '/chromosomes.fasta'
    pla_references_fasta = eval_dir + '/plasmids.fasta'
    subprocess.call('rm -f %s %s' % (chr_references_fasta, pla_references_fasta), shell = True)
    for acc in chromosomes:
        file = '%s/%s.fasta' % (out_dir, acc)
        subprocess.call('bash %s %s %s; cat %s >> %s; rm %s' % (extract_script, acc, file, file, chr_references_fasta, file), shell = True)
    for acc in plasmids:
        file = '%s/%s.fasta' % (out_dir, acc)
        subprocess.call('bash %s %s %s; cat %s >> %s; rm %s' % (extract_script, acc, file, file, pla_references_fasta, file), shell = True)

    blast_db = eval_dir + '/plasmid_references'
    subprocess.call('makeblastdb -in %s -dbtype nucl -out %s; ' % (pla_references_fasta, blast_db), shell = True)

    #blast_db = eval_dir + '/chr_references'
    #subprocess.call('makeblastdb -in %s -dbtype nucl -out %s; ' % (chr_references_fasta, blast_db), shell = True)    
    
    quast_labels = []
    quast_files = []
    
    
    ## evaluate MILP-based plasmids
    # first, map plasmid genes to Unicycler assembly and filter the mapping
    print('Performing evaluations with MILP strategy...')
    subprocess.call('mkdir %s/MILP' % eval_dir, shell = True)
    
    # MILP / partial training                     
    MILP_putative_queries = out_dir + '/MILP/putative_plasmids.fasta'
    MILP_putative_mapping = eval_dir + '/MILP/MILP_putative_map.csv'
    #greedy_questionable_queries = out_dir + '/greedy/plasmids/greedy/questionable_plasmids.fasta'
    #greedy_questionable_mapping = eval_dir + '/greedy/greedy_questionable_map.csv'
    subprocess.call( 'blastn -task megablast -query %s -db %s -out %s -outfmt 6; ' % (MILP_putative_queries, blast_db, MILP_putative_mapping) \
                    #+ 'blastn -task megablast -query %s -db %s -out %s -outfmt 6; ' % (greedy_questionable_queries, blast_db, greedy_questionable_mapping) \
                    + 'python %s %s %s %s %s/MILP/MILP_eval.csv -i "putative;%s;%s" >> %s; ' % (compare_script, sample_id, chr_references_fasta, pla_references_fasta, eval_dir, MILP_putative_queries, MILP_putative_mapping, results) , shell = True)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument("proj_dir", help='path to PlasEval dir')
    parser.add_argument('sample', help = 'Name/ID of sample')
    parser.add_argument('out_dir', help = 'path to output directory')
    parser.add_argument('results', help = 'CSV file of scores')
    args = parser.parse_args()

    proj_dir = args.proj_dir
    compare_script = proj_dir + '/eval_scripts/compare_plasmids.py'
    extract_script = proj_dir + '/data/2018-05-17__MOB-suite_benchmark/extract.sh'
    sample_script = proj_dir + '/data/2018-05-23__MOB-suite_benchmark_reads/sample.sh'
    
    evaluate(args.sample, args.out_dir, args.results, compare_script, extract_script, sample_script)
    print('Finished analysis of sample %s.' % args.sample)

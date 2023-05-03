# PlasEval
## Comparison script
### Input
1. Two .tsv files with predicted / reference plasmids as sets of contigs. The file should contain one chain per line.<br/>
Format:<br/>
Plasmid_ID	Contig_ID 	Contig_Length<br/>
2. Path to output directory and name of output file.

#### Usage
```
python plasmid_comparison_main.py --l left_plasmids.tsv --r right_plasmids.tsv --out out_dir --res out_file
```

### Output
A log file that contains the following information for each possible labelling combination. The combinations are listed in increasing order of the dissimilarity score.:<br\>
1. Combination number,<br/>
2. Dissimilarity score,<br/>
3. Total length of contigs that are unique to one set of plasmids and half the total length of contigs that are common but grouped differently in both sets of plasmids,<br/>
4. Labelling of contigs in both sets of plasmids,<br/>
5. Splits in both sets of plasmids

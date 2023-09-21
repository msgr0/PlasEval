# PlasEval
## Comparison script
### Input
1. Two TSV files with predicted / reference plasmids as sets of contigs. The file should contain one chain per line.<br/>
Format:<br/>
Plasmid_ID	Contig_ID 	Contig_Length<br/>
2. Path to output file and log file.

#### Usage
```
python plasmid_comparison_main.py --l LEFT_BINS_TSV --r RIGHT_BINS_TSV --out_file OUT_FILE --log_file LOG_FILE
```
where `l` and `r` are TSV files, each with one set of plasmid bins. `out_file` is the path to the output file while `log_file` is the path to the log file.

### Output
The output file contains the following information:<br\>
1. Cumulative length of contigs present in at least one of set of plasmid bins,<br/>
2. Cost of cuts: splitting bins from first set of plasmid bins,<br/>
3. Cost of joins: splitting bins from second set of plasmid bins,<br/>
4. Cumulative length of contigs present only in the first set,<br/>
5. Cumulative length of contigs present only in the second set,<br/>
6. Dissimilarity score
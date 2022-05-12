# PlasEval

## Input
1. A file with predicted contig chains with or without orientations. The file should contain one chain per line.<br/>
Format with orientation:<br/>
plasmid_1;27+,28-,12+<br/>
plasmid_2;13-,1+,4+<br/>

Format without orientation:<br/>
plasmid_1;27,28,12<br/>
plasmid_2;13,1,4<br/>

2. A file that maps a contig to a plasmid, chromosome or ambiguous sequence with one contig per line. <br/>
Format:<br/>
short_read_contig_id,short_read_contig_length,mapped_against_hybrid_contig_id,hybrid_contig_length,number_of_residue_matches,is_circular,label

3. Path to output directory and name of output file.

### Usage
```
python evaluate_sample.py --pred prediction_file --map mapping_file --out out_dir --res out_file --amb 0/1 --ori 0/1
```

Additional arguments
```
--amb			Argument to indicate if contigs mapped ambiguous sequences should be considered plasmidic. If yes, then pass value 1, else 0.
--ori			Argument to indicate if predicted contig chains have oriented contigs. If yes, then pass value 1, else 0.                           
```
# PlasEval

**TO DO**: quick intro

## Plasmid binning

**TO DO [cc]**: quick description of plasmid binning

## Plasmid binning evaluation

PlasEval is a tool designed for evaluating the results of plasmid binning methods. It has two modes: evaluate and compare. 
1. Evaluate: This mode takes a set of predicted plasmid bins as well as a set of ground truth plasmid bins as input. It then computes the precision and recall statistics for each predicted plasmid bin and ground truth bin respectively. <br/>
2. Compare: This mode compares two sets of plasmid bins. It computes a dissimilarity score between the two sets of plasmid bins based on transforming one of the sets into the other using a sequence of splits and joins, weighted by contigs length.

**TO DO [cc]**: quick description of the dissimilarity.

## Usage

### Input: plasmid bins file

The main input in both modes are the two plasmid bin files. These are two TSV files with predicted / ground truth plasmids as sets of contigs. The file should have a header row with the column names, specifically the plasmid name, contig id and contig length, as shown in the example below. The file should contain one contig per line with contig information provided under the respective columns. The file can contain other information as long as the three columns mentioned above are provided.<br/>
plasmid	contig 	contig_len<br/>
P1	C1 	2000<br/>
P1	C2 	3000<br/>
P2	C1 	2000<br/>
P2	C3	4000<br/>
P2	C1	2000<br/>
...<br/>

### Input: numeric parameters

Length threshold: Contigs below the length threshold are excluded from the computation of precision and recall. This is an optional input with a default value of 0. 
Maximum number of calls.

### Input: output files

#### Usage
```
python eval/evaluate_bins.py --pred PREDICTED_BINS_TSV --gt GROUNDTRUTH_BINS_TSV --out OUT_DIR --res OUT_FILE --thlen LEN_THRESHOLD
```
where `pred` and `gt` are TSV files, with the set of predicted and ground truth plasmid bins respecitvely. `out` is the path to the output directory while `res` is the name of the output file. `thlen` is the integer length threshold as mentioned above.

3. Compare mode: <br/>
In addition to the two plasmid bin files, the compare mode requires the path to output file and log file as input.

#### Usage
```
python compare/plasmid_comparison_main.py --l LEFT_BINS_TSV --r RIGHT_BINS_TSV --out_file OUT_FILE --log_file LOG_FILE
```
where `l` and `r` are TSV files, each with one set of plasmid bins. `out_file` is the path to the output file while `log_file` is the path to the log file.

### Output
1. The output file for the evaluate mode contains the following information:<br/>
	a. Number and names of predicted plasmid bins <br/>
	b. Precision details: For each predicted bin, for both contig level and basepair level precision, the name of the best matched ground truth bin and the corresponding precision values <br/> 
	c. Recall details: For each ground truth bin, for both contig level and basepair level recall, the name of the best matched predicted plasmid bin and the corresponding recall values <br/> 
	d. Overall contig level and basepair level statistics

2. The output file for the compare mode contains the following information:<br/>
	a. Cumulative length of contigs present in at least one of set of plasmid bins,<br/>
	b. Cost of cuts: splitting bins from first set of plasmid bins,<br/>
	c. Cost of joins: splitting bins from second set of plasmid bins,<br/>
	d. Cumulative length of contigs present only in the first set,<br/>
	e. Cumulative length of contigs present only in the second set,<br/>
	f. Dissimilarity score

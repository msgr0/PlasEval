# PlasEval

PlasEval is a tool aimed at evaluating the accuracy and at comparing methods for the problem of **plasmid binning**.

## Plasmid binning

Plasmid binning aims at detecting, from the draft assembly of a single isolate bacterial pathogen, groups of contigs assumed eachto originate from a plasmid present in the sequenced isolate. Recent methods for plasmid binning include <a href="https://github.com/cchauve/PlasBin-flow">PlasBin-Flow</a>, <a href="https://github.com/phac-nml/mob-suite">MOB-recon</a> and <a href="https://gitlab.com/sirarredondo/gplas">gplas</a>.

The result of applying a plasmid binning method to a draft assembly is a set of unordered groups of contigs, called **plasmid bins**.

## Plasmid binning evaluation and comparison

### Evaluation
In its **evaluation** mode of PlasEval, given a set of *predicted plasmid bins* resulting from a plasmid binning tools and a *ground truth* set of plasmid bins, where each true plasmid bin contains all the contigs that belong to one of the plasmid present in the sequenced isolate, PlasEval computes three statistics, the *precision*, the *recall* and the *F1-score* of the predicted plasmid bins with respect to the ground truth.

For a group $X$ of contigs, we denote by $L(X)$ the cumulated length of the contigs in $X$. 
For a predicted plasmid bin $P$ and a ground truth plasmid bin $T$, we define the *overlap* between $P$ and $T$, denoted by $o(P,T)$, as the cumulated length of the contigs presents in both $P$ and $T$.
Given a set $\cal{P}$ of predicted plasmid bins and a set $\cal{T}$ of ground truth plasmid bins, w define the precision $p(\cal{P},\cal{T})$ and recall $r(\cal{P},\cal{T})$ as follows:
$$p(\cal{P},\cal{T}) = \frac{\sum\limits_{P\in\cal{P}} \max\limits_{T\in\cal{T}} o(P,T)}{\sum\limits_{P\in\cal{P}}L(P)}$$
$$r(\cal{P},\cal{T}) = \frac{\sum\limits_{T\in\cal{T}} \max\limits_{P\in\cal{P}} o(P,T)}{\sum\limits_{T\in\cal{T}}L(T)}$$
and the F1-score as the arithmetic mean of the precision and recall
$$F_1(\cal{P},\cal{T}) = 2\frac{{p}(\cal{P},\cal{T}){r}(\cal{P},\cal{T})}{{p}(\cal{P},\cal{T})+{r}(\cal{P},\cal{T})}.$$

### Comparison
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

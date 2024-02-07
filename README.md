# PlasEval

PlasEval is a tool aimed at evaluating the accuracy and at comparing methods for the problem of **plasmid binning**.

## Plasmid binning

Plasmid binning aims at detecting, from the draft assembly of a single isolate bacterial pathogen, groups of contigs assumed eachto originate from a plasmid present in the sequenced isolate. Recent methods for plasmid binning include <a href="https://github.com/cchauve/PlasBin-flow">PlasBin-Flow</a>, <a href="https://github.com/phac-nml/mob-suite">MOB-recon</a> and <a href="https://gitlab.com/sirarredondo/gplas">gplas</a>.

The result of applying a plasmid binning method to a draft assembly is a set of unordered groups of contigs, called **plasmid bins**.

## Plasmid binning evaluation and comparison

### Evaluation
In the **evaluation** mode of PlasEval, given a set of *predicted plasmid bins* resulting from a plasmid binning tools and a *ground truth* set of plasmid bins, where each true plasmid bin contains all the contigs that belong to one of the plasmid present in the sequenced isolate, PlasEval computes three statistics, the *precision*, the *recall* and the *F1-score* of the predicted plasmid bins with respect to the ground truth.

For a group $`X`$ of contigs, we denote by $`L(X)`$ the cumulated length of the contigs in $`X`$. 
For a predicted plasmid bin $`P`$ and a ground truth plasmid bin $`T`$, we define the *overlap* between $`P`$ and $`T`$, denoted by $`o(P,T)`$, as the cumulated length of the contigs presents in both $`P`$ and $`T`$.
Given a set $`A`$ of predicted plasmid bins and a set $`B`$ of ground truth plasmid bins, we define the precision $`p(A,B)`$ and recall $`r(A,B)`$ as follows:
```math
p(A,B) = \frac{\sum\limits_{P\in A} \max\limits_{T\in B} o(P,T)}{\sum\limits_{P\in A}L(P)}, r(A,B) = \frac{\sum\limits_{T\in B} \max\limits_{P\in A} o(P,T)}{\sum\limits_{T\in B}L(T)}.
```  

The F1-score is the arithmetic mean of the precision and recall
```math
F_1(A,B) = 2\frac{{p}(A,B){r}(A,B)}{{p}(A,B)+{r}(A,B)}.
```  

### Comparison
In the **comparison** mode of PlasEval, given two sets of *predicted plasmid bins* (either resulting from two plasmid binning tools or from a plasmid binning tool and a ground truth), PlasEval computes three statistics, PlasEval computes a *dissimilarity measure* that indictaes how much both sets of predicted plasmid bins are in agreement. 
The full details of the dissimilarity measure computed by PlasEval are available in the preprint <a href="">PlasEval: a framework for comparing and evaluating plasmid binning tools</a> and we describe below its general principle.
The dissimilariy value is the sum of 4 terms accounting respectively for
- *extra contigs*: the contigs present in $`A`$ but not in $`B`$;
- *missing contigs*: the contigs present in $`B`$ but not in $`A`$;
- *splits*: splitting of the plasmid bins of $`A`$ to obtain a set of intermediate bins defined as the intersection of $`A`$ and $`B`$;
- *joins*: joining the splitted plasmid bins into the plasmid bins of $`B`$.

Each of these components incur a *cost*, parameterized by a parameter $`\alpha \in [0,1]`$:
- each extra or missing contig $`c`$ results in a cost $`\ell(c)^\alpha`$, where $`\ell(c)`$ is the length of $`c`$;
- each split of a plasmid bin $`P`$ into two smaller bins $`P',P''`$ results in a cost $`\min(L(P')^\alpha,L(P'')^\alpha)`$;
- each join of two intermediate plasmid bins $`P',P''`$ into  larger bin $`P`$ results in a cost $`\min(L(P')^\alpha,L(P'')^\alpha)`$.
- 
So with $`\alpha=0`$, the length of contigs is not accounted for and only set-theoretic operations define the dissimilarity value, while with $`\alpha=1`$ it is fully accounted for.
By default $`\alpha=0.5`$.

In the case where some contigs are *repeated*, i.e. a contig appears in more than one plasmid bin of $`A`$ and/or $`B`$, a *branch-and-bound* algorithm computes the pairing between repeats contigs in $`A`$ and in $`B`$ that results in the minimum dissimilarity value.

Finally the dissimilarity obtained as described above is *normalized* into a value in $`[0,1]`$ by dividing it by 
```math
\sum\limits_{P\in A}\sum\limits_{c\in P} \ell(c)^{\alpha} + \sum\limits_{Q\in B}\sum\limits_{c\in Q} \ell(c)^{\alpha}.
```
  
## Usage

### Input: plasmid bins file

The main input in both modes ofPlasEval are the two plasmid bins files. 
A plasmid bins file is a TSV file that describes a set of plasmids bins in a forma where each row contains three pieces of information: the identifier of a plasmid bin, the identifier of a contig that belongs to this plasmid bin and the length of the contig.
The file should have a header row with column names `plasmid`, `contig`, `contig_len`. 
An example is provided below that describes the set of plasmid bins where bin `P1` contains contigs `C1,C2` and bin `P2` contains contigs `C1`, `C3` and `C4`.
```
plasmid	contig 	contig_len
P1	C1 	2000
P2	C3 	3000
P1	C2 	2000
P2	C1	2000
P2	C4	2000
```

If a contig appears in several copies in a plasmid bin, only one copy is accounted for.
Other information can be included in the TSV file, but will not be accounted for by PlasEval.

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

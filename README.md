# [software] TRIBECaller

![img](https://yy4h2ftvat.feishu.cn/space/api/box/stream/download/asynccode/?code=NTFkMTBiNDcwNTAxMTAwMGI0NzVlNDRlYWQwNDUzODFfMGg1SG12cVJVZWc5VjB1V1NEcTUyYWFIdHU2djROakVfVG9rZW46Ym94Y25qenprWktTOWZkZG55bG9DdUpGU2plXzE2MTE1NTY2MDE6MTYxMTU2MDIwMV9WNA)

## **Project Timeline**

Unable to paste block outside Docs

## **Project Description**

#### **Background**

RNA editing are ubiquitous in a wide range of organisms and it is an crucial post-transcriptional mechanism to regulate the function of primary mRNA through insertion, deletion, or modification (editing) of specific nucleotides. RNA editing have long been known to occur in tRNAs, rRNAs, and mRNAs. Two common types of RNA editing involve deamination reaction, either by deamination of cytidine (C) produces uridine (U), or deamination of adenine (A) to inosine (I). In mammals, two adenosine deaminase (ADAR) proteins have been found to catalyzes adenine-to inosine (A-to-I) conversion in dsRNA, without additional factors. The inosine can be then converted to guanosine and paired with cytosine. 

![img](https://yy4h2ftvat.feishu.cn/space/api/box/stream/download/asynccode/?code=YThiYjZhYTJkNzE2ZjA0YjYxYTQwNzg3Mjc0ZmM4MzhfUzdKQndBWXgzQWFjYlhZenJqdlh3ckllc2pWSU83NDJfVG9rZW46Ym94Y24yWUlNN09hVm8zejVoVVlQc1I4SnlmXzE2MTE1NTY2MDE6MTYxMTU2MDIwMV9WNA)

(Figure derived from RNA (2001) Dieter Söll, Susumu Nishimura and Peter Moore (Eds.)

*Drosophila melanogaster* has a single Adar gene encoding a protein related to mammalian ADAR2 that edits adenine in early mRNA transcripts. Similar to the mammalian ADAR2, endogenous Drosophila ADAR is a modular enzyme consisting a dsRNA binding motif and a catalytic domain. McMahon et al. (2016) have developed TRIBE (**T**argets of **R**NA-binding proteins **I**dentified **B**y **E**diting) using only the fusion protein containing RNA binding protein of interest and the catalytic domain of *Drosophila* ADAR (ADARcd). The irreservable editing would allow RNA-sequencing to identify the editing sites and thus the binding targets of RBPs will be discovered.

![img](https://yy4h2ftvat.feishu.cn/space/api/box/stream/download/asynccode/?code=ZjNlNDc2NGQ1MGFmNzg3ZWM5MDJjNGM1ZmIyMjkzN2JfNWhiVTVYMTNVWEVqTkpoeGdmZHlaNGtLMGRucmh5RlZfVG9rZW46Ym94Y253SXVIRDBiZUhjTXo1NGpiNWR2bDJiXzE2MTE1NTY2MDE6MTYxMTU2MDIwMV9WNA)

(Figure derived from McMahon et al. (2016))

Scripts or software written in different programming languages have provided tools to identify these RNA editing sites, but the performance and scalability of these methods require further improvement. 

### Rationale

**TRIBECaller** is a software written in Python to find editing sites using raw data (fastq) from RNA-sequencing data of TRIBE experiment.

**TRIBECaller** converts every read from the bam file to the forward strand of the reference genome, and restricts A-to-G mutations in transcripts encoded by the forward strand and T-to-C mutations in transcripts encoded by the reverse strand.

![img](https://yy4h2ftvat.feishu.cn/space/api/box/stream/download/asynccode/?code=MDBhY2RlYmUzNzQwZDliYzNmMTVmMjZmYmJmMWQxOWFfcEhscTV2cHRuTnpzQ3lta2liS0xsTzdCbExXazVmeUdfVG9rZW46Ym94Y25hR2NtT0FBeWllNUhnWUJiS0RzQWVlXzE2MTE1NTY2MDE6MTYxMTU2MDIwMV9WNA)

Currently, **TRIBECaller** use fisher's exact test to call RNA editing sites.

|  | **A(forward) or T(reverse)** | **G(forward) or C(reverse)** | **Row Total** |
| ---------------- | ---------------------------- | ---------------------------- | ------------- |
| **Experiment**   | m| n| m+n   |
| **Control**  | a| b| a+b   |
| **Column Total** | m+a  | n+b  | m+n+a+b   |

Fisher showed that the probability of obtaining any such set of values was given by the hypergeometric distribution:

$$p = \frac{(m+n)!(a+b)!(m+a)!(n+b)!}{m!n!a!b!(m+n+a+b)!}$$, $$\mathrm{odds} = \frac{\frac{n}{m}}{\frac{b}{a}}$$

The program will also output the value of difference of G/C content defined by:

$$\mathrm{diff} = \frac{n}{m+n} - \frac{b}{a+b}$$

### Implementation

**TRIBECaller** is implemented in Python3, dependent on a python wrapper of samtools, pysam. Parallel computation of editing sites is supported.

Using mapReduce as our parallel computation model, **TRIBECaller** uses multiple worker threads to compute the editing sites and test statistics.

![img](https://yy4h2ftvat.feishu.cn/space/api/box/stream/download/asynccode/?code=ZTM1ODM4YmU3ODU5NmEyMGVlOTZkNTYyZDcyM2ExYWZfdTBEVmM1U1dwamx1ejFQbnJ5M2d2cDI4RVBSUW5mRlpfVG9rZW46Ym94Y25HcHhtSjRRMEI4YzdlbUxzWTFtRmVoXzE2MTE1NTY2MDE6MTYxMTU2MDIwMV9WNA)

#### **Performance**

In one test, calling editing events with **TRIBECaller** using an 1.2G experiment bam file and 961M control bam file only costs 120 minutes, using 12 parallel threads, Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz.

## **User Document**

**TRIBECaller** is currently under development and testing. You can try to use version 0.0.1 on zje610 server: /home/zje610/workspace/labW/Snowxue/TRIBECaller or wanlu's hoffman: ~/wanluliu/TRIBE/TRIBECaller

The software is temporarily called by python main.py. Binary executable will be available in later versions. **TRIBECaller** contains major calling modules and two plotting module to visualize the result.

```
python main.py  -v

TRIBECaller 0.0.1



python main.py -h

usage: main.py [-h] [-v]

   {callEditingSites,plotEditingSite,plotEditingRegion,computeCoverage}

   ...



Editing site caller for TRIBE



positional arguments:

  {callEditingSites,plotEditingSite,plotEditingRegion,computeCoverage}

callEditingSitesMain TRIBECaller function

plotEditingSite Plot nucleotides from reads within a genomic region

plotEditingRegion   Plot editing events within a genomic region

computeCoverage Compute coverage and nucleotides content



optional arguments:

  -h, --helpshow this help message and exit

  -v, --version show program's version number and exit



For command line options for each command, type COMMAND -h
```

## Subprograms

### callEditingSites

```
python main.py plotEditingSite -husage: main.py callEditingSites [-h] -t TARGET -c CONTROL [-b BINSIZE] -o

OUTPREFIX [-contig CONTIG] [-gz] [-p THREAD]



optional arguments:

  -h, --helpshow this help message and exit

  -t TARGET, --target TARGET

input target bam file

  -c CONTROL, --control CONTROL

input control bam file

  -b BINSIZE, --binSize BINSIZE

bin size of calling editing event

  -o OUTPREFIX, --outPrefix OUTPREFIX

prefix of output bed file

  -contig CONTIG, --contig CONTIG

chromosome of interests

  -gz, --gzip   output gzipped file

  -p THREAD, --thread THREAD

number of threads to use



Examples: python main.py callEditingSites -t Experiment_rep1_1_srt.bam -c

Control_rep1_1_srt.bam -o testPrefix
```

callEditingSites calls all editing sites according to an input experiment (TRIBE) bam file and a control experiment bam file. The output will be a bed file containing all coordinates of the editing events.

**Arguments**

-t/--target: Input target bam file

-c/--control: Input control bam file

-o/--out: the prefix of output (ended with .bed extension)

-contig: if provided, only the reads mapped within these contig/chromosome will be calculated. If multiple contig are provided, they should be comma separated.

-gz/--gzip: if setted, the program will output to a gzipped file.

-p/--thread: number of threads to use in the program.

**Example output**

```
1   20249   20250   0.4166666666666667  6   6   11  1   0.034324942791762014-

1   134346  134347  0.752   6   4   0   0.030303030303030304-

1   135258  135259  0.05405405405405406 35  2   160 0   0.03449704754998473 +

1   136217  136218  0.8333333333333334  1   5   12  0   0.0007002801120448174   -

1   136997  136998  0.6 2   3   8   0   0.03496503496503495 -

1   137634  137635  0.05292297671389024 83  6   136 2   0.04220621740014615 -

1   144940  144941  0.6666666666666666  1   2   12  0   0.028571428571428577+

1   149037  149038  1.0 0   7   2   0   0.027777777777777794-

1   190358  190359  0.4 3   2   16  0   0.0476190476190477  +

1   195419  195420  1.0 0   3   6   0   0.011904761904761908-
```

**Column 1:** Chromosome name

**Column 2:** start of the position

**Column 3:** end of the position

**Column 4:** Difference of nucleotides composition. If the editing event comes from the forward strand, the value will be the difference of (G/(A+G) between the experiment and control RNA seq. If the editing event comes from the reverse strand, the value will be the difference of (C/(T+C) between the experiment and control RNA seq.

**Column 5:** Count of A/T in the experiment sample

**Column 6:** Count of G/C in the experiment sample

**Column 7:** Count of A/T in the control sample

**Column 8:** Count of G/C in the controlsample

**Column 9:** p-value of the fisher-exact test

**Column 10** : Strand. + for the forward strand and - for the reverse strand.

#### computeCoverage

```
python main.py computeCoverage -h

usage: main.py computeCoverage [-h] -t TARGET -c CONTROL -o OUTPREFIX

   [-contig CONTIG] [-gz] [-p THREAD]



optional arguments:

  -h, --helpshow this help message and exit

  -t TARGET, --target TARGET

input target bam file

  -c CONTROL, --control CONTROL

input control bam file

  -o OUTPREFIX, --outPrefix OUTPREFIX

prefix of output bed file

  -contig CONTIG, --contig CONTIG

chromosome of interests

  -gz, --gzip   output gzipped file

  -p THREAD, --thread THREAD

number of threads to use



Examples: python main.py computeCoverage -t Experiment_rep1_1_srt.bam -c

Control_rep1_1_srt.bam -o testPrefix
```

computeCoverage compute nucleotide (A,T,C,G) composition from every nucleotide position. It will output a .txt file, either uncompressed or compressed by gzip. 

-t/--target: Input target bam file

-c/--control: Input control bam file

-o/--out: the prefix of output (ended with .bed extension)

-contig: if provided, only the reads mapped within these contig/chromosome will be calculated. If multiple contig are provided, they should be comma separated.

-gz/--gzip: if setted, the program will output to a gzipped file.

-p/--thread: number of threads to use in the program.

**Example output:**

```
1   14404   14405   0   2   0   0   2   0   1   0   0   1

1   14405   14406   0   2   0   0   2   0   1   0   0   1

1   14406   14407   0   0   2   0   2   0   0   1   0   1

1   14407   14408   0   2   0   0   2   0   1   0   0   1

1   14408   14409   0   0   0   2   2   0   0   0   2   2

1   14409   14410   0   0   6   0   6   0   0   2   0   2

1   14410   14411   0   6   0   0   6   0   2   0   0   2

1   14411   14412   0   0   6   0   6   0   0   2   0   2

1   14412   14413   6   0   0   0   6   2   0   0   0   2

1   14413   14414   0   0   0   7   7   0   0   0   3   3
```

**Column 1:** Chromosome name

**Column 2:** start of the position

**Column 3:** end of the position

**Column 4-7:** count of A,T,C,G in the experiment sample

**Column 8:** total count of reads in this position in the experiment sample

**Column 9-12:** count of A,T,C,G in the control sample

**Column 13:**  total count of reads in this position in the control sample

#### plotEditingSite

After calling the editing events, you may want to visualize the nucleotide change in specific region undergoing editing. plotEditingSite plots the nucleotides percentage in a specific site.

```
python main.py plotEditingSite -h

usage: main.py plotEditingSite [-h] -t TARGET -c CONTROL -r REGION -o OUT

   [--dpi DPI]



optional arguments:

  -h, --helpshow this help message and exit

  -t TARGET, --target TARGET

input taregt bam file

  -c CONTROL, --control CONTROL

input control bam file

  -r REGION, --region REGION

genomic region of interests. Format: 1:1072894-1078023

  -o OUT, --out OUT prefix of output figure

  --dpi DPI dpi of the output figure



Examples: python main.py plotEditingSite -t Experiment_rep1_1_srt.bam -c

Control_rep1_1_srt.bam --region 1:1214232:1214243
```

### plotEditingRegion

After calling the editing events, you may want to visualize the reads coverage and number of editing events in a genomic region. You can visualize the reads and editing events in the range of a gene.

```
python main.py plotEditingRegion -h

usage: main.py plotEditingRegion [-h] -t TARGET -c CONTROL [-r REGION]

 [-g GENE] -gtf GTF -o OUT [--dpi DPI]



optional arguments:

  -h, --helpshow this help message and exit

  -t TARGET, --target TARGET

input taregt bam file

  -c CONTROL, --control CONTROL

input control bam file

  -r REGION, --region REGION

genomic region of interests. Format: 1:1072894-1078023

  -g GENE, --gene GENE  gene of intersets. For multiple genes, use comma

separated input.

  -gtf GTF, --gtf GTF   gtf file to plot the genes.

  -o OUT, --out OUT prefix of output figure

  --dpi DPI dpi of the output figure

  
```

![img](https://yy4h2ftvat.feishu.cn/space/api/box/stream/download/asynccode/?code=ZmY2YTk4YTJjMTRhNGM3OTk3MWQxYjIzMGNjMWU0YWFfaUE5U0cwNUNBMU5CWEx0dE8wc3ZwNUdqUmdZeU5LRDZfVG9rZW46Ym94Y25OSVF4QTRKeEpwc2YzM1FGMHVaSE5mXzE2MTE1NTY2MDE6MTYxMTU2MDIwMV9WNA)

For example,

```
python main.py plotEditingRegion -t /Users/snowxue/Documents/TribeTest/CCT4_rep1_1_CCT4_rep1_2Aligned.sortedByCoord.sf4_srt.bam -c /Users/snowxue/Documents/TribeTest/H9_rep1_1.f_H9_rep1_2.fAligned.sortedByCoord.sf4_srt.bam --gene DFFA -gtf ~/Documents/refData/Homo_sapiens.GRCh38.97.chr.gtf -o test.pdf
```

gives

![img](https://yy4h2ftvat.feishu.cn/space/api/box/stream/download/asynccode/?code=NGNhNGMwOTkwZTVjZDRlZTNmODg1MzhiODc2ZWNhNDZfRXpPbTdGc1BSdlpqZ0hyUjJZRlRvbk1PQ05HaXAxa0dfVG9rZW46Ym94Y25TVW1KaUJ6VFY1clBIWEVYTDFEQzJmXzE2MTE1NTY2MDE6MTYxMTU2MDIwMV9WNA)

The first row is the gene annotation

The second row is the read coverage in the experiment sample

The third row is the computed editing sites

The fourth row shows the "diff" values

Tfirth row is the read coverage in the control sample

# **Release Note**

#### **0.0.0** 

[**2021.1.18]** Initial release.

#### **0.0.1**

[**2021.1.24]** Multithreading and MapReducing to achieve faster computation.

### **Future Development Plan**

- Parallel computation of the A(G)T(C) dictionary (Advanced data structure might be used; Try to overcome the GIL)
- MapReduce mechanism on computing the A(G)T(C) dictionary might be considered.
- Visualization of editing events in a wider genomic range.
- Fetching the gtf file is slow. However, samtools do not support indexing the gtf format. We should make a custom index in the later versions.
- Deal with stranded library.
- Compile binary executables.
- Whether to use beta-binomial distribution described in Nguyen et al., 2020.
- High memory mode for faster computation.
- More command line arguments, for customized editing site calling.
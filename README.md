# [software] TRIBECaller

## Project Description

TRIBECaller is a software written in Python to find editing sites using data from TRIBE experiment.

![img](https://yy4h2ftvat.feishu.cn/space/api/box/stream/download/asynccode/?code=NDM3NTYwNDk0YTFiYzhlY2Q2Y2Y0MWQxNzRkNzkwNmZfUUgxd2hyYlJ0aDNqcnhadEM2NUZ3NDRtZjZLV3hJZ3ZfVG9rZW46Ym94Y25vNXBmemEwd0pUVlF6bGYwNjRtY09FXzE2MTA5Mzc4ODc6MTYxMDk0MTQ4N19WNA)

## User Document

TRIBECaller is currently under development and testing. You can try to use the demo version (0.0.0) on zje610 server: /home/zje610/workspace/labW/Snowxue/TRIBECaller

The software is temporarily called by python main.py. Binary executable will be available in later versions.

```
$ python main.py  -v

TRIBECaller 0.0.0

$ python main.py  -h

usage: main.py [-h] [-v]

               {callEditingSites,plotEditingSite,plotEditingRegion} ...



Editing site caller for TRIBE



positional arguments:

  {callEditingSites,plotEditingSite,plotEditingRegion}

    callEditingSites    Main TRIBECaller function

    plotEditingSite     Plot nucleotides from reads within a genomic region

    plotEditingRegion   Plot editing events within a genomic region



optional arguments:

  -h, --help            show this help message and exit

  -v, --version         show program's version number and exit



For command line options for each command, type COMMAND -h
```

### Subprograms

#### callEditingSites

```
$ python main.py plotEditingSite -h

usage: main.py plotEditingSite [-h] -t TARGET -c CONTROL -r REGION -o OUT

                               [--dpi DPI]



optional arguments:

  -h, --help            show this help message and exit

  -t TARGET, --target TARGET

                        input target bam file

  -c CONTROL, --control CONTROL

                        input control bam file

  -r REGION, --region REGION

                        genomic region of interests. Format: 1:1072894-1078023

  -o OUT, --out OUT     prefix of output figure

  --dpi DPI             dpi of the output figure



Examples: python main.py plotEditingSite -t Experiment_rep1_1_srt.bam -c

Control_rep1_1_srt.bam --region 1:1214232:1214243
```

callEditingSites calls all editing sites according to an input experiment (TRIBE) bam file and a control experiment bam file. The output will be a bed file containing all coordinates of the editing events.

**Arguments**

-t/--target: Input target bam file

-c/--control: input control bam file

-r/--region: Under development

-o/--out: prefix of output figure in pdf

#### plotEditingSite

After calling the editing events, you may want to visualize the nucleotide change in specific region undergoing editing. plotEditingSite plots the nucleotides percentage in a specific site.

```
$ python main.py plotEditingSite -h

usage: main.py plotEditingSite [-h] -t TARGET -c CONTROL -r REGION -o OUT

                               [--dpi DPI]



optional arguments:

  -h, --help            show this help message and exit

  -t TARGET, --target TARGET

                        input taregt bam file

  -c CONTROL, --control CONTROL

                        input control bam file

  -r REGION, --region REGION

                        genomic region of interests. Format: 1:1072894-1078023

  -o OUT, --out OUT     prefix of output figure

  --dpi DPI             dpi of the output figure



Examples: python main.py plotEditingSite -t Experiment_rep1_1_srt.bam -c

Control_rep1_1_srt.bam --region 1:1214232:1214243
```

#### plotEditingRegion

After calling the editing events, you may want to visualize the reads coverage and number of editing events in a genomic region. You can visualize the reads and editing events in the range of a gene.

```
$ python main.py plotEditingRegion -h

usage: main.py plotEditingRegion [-h] -t TARGET -c CONTROL [-r REGION]

                                 [-g GENE] -gtf GTF -o OUT [--dpi DPI]



optional arguments:

  -h, --help            show this help message and exit

  -t TARGET, --target TARGET

                        input taregt bam file

  -c CONTROL, --control CONTROL

                        input control bam file

  -r REGION, --region REGION

                        genomic region of interests. Format: 1:1072894-1078023

  -g GENE, --gene GENE  gene of intersets. For multiple genes, use comma

                        separated input.

  -gtf GTF, --gtf GTF   gtf file to plot the genes.

  -o OUT, --out OUT     prefix of output figure

  --dpi DPI             dpi of the output figure

  
```

![img](https://yy4h2ftvat.feishu.cn/space/api/box/stream/download/asynccode/?code=ZWRkYTQyYWEwNjU1ZGQxN2RiMjkxYWQzNDM2ZGNjNjNfOU1XZ2g2OVhnV3ZMekkzREFaN1VTUGlTWlkya1Y3QnVfVG9rZW46Ym94Y25OSVF4QTRKeEpwc2YzM1FGMHVaSE5mXzE2MTA5Mzc4ODc6MTYxMDk0MTQ4N19WNA)

For example,

```
$ python main.py plotEditingRegion -t /Users/snowxue/Documents/TribeTest/CCT4_rep1_1_CCT4_rep1_2Aligned.sortedByCoord.sf4_srt.bam -c /Users/snowxue/Documents/TribeTest/H9_rep1_1.f_H9_rep1_2.fAligned.sortedByCoord.sf4_srt.bam --gene DFFA -gtf ~/Documents/refData/Homo_sapiens.GRCh38.97.chr.gtf -o test.pdf
```

gives

![img](https://yy4h2ftvat.feishu.cn/space/api/box/stream/download/asynccode/?code=YzJiMzRhMDI0NTE5NzVkNTNlMWJlYWY0NDBjYzYzYWRfbUJZZk1nZ1NPU2F4Z01lMlJjc3dCN3o5M3RNUGtZQ2tfVG9rZW46Ym94Y25kZVdzV1NBdFYzTnA1Q1ZyVEpHZ3FjXzE2MTA5Mzc4ODc6MTYxMDk0MTQ4N19WNA)

## Release Note

- 0.0.0 2021.1.18 Initial release
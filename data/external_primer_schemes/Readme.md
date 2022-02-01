# source
* V1200 bp scheme from: https://zenodo.org/record/3897530#.X4QxOpqxUay
* Primers V to V4* -> https://github.com/artic-network/artic-ncov2019

They are just copied here for workflow stability and ease of use.

# IMPORTANT
+ PRIMER DIRS NEED TO START WITH A "V"
+ BED FILES structure
```
MN908947.3	22	46	varskip-0317-1_01_LEFT	nCoV-2019_1	+
```
* first column MN908947.3 (the fasta header of the reference)
* 4th colum sort it (sort -k4)
    * also remove _alt   might be problematic
* 5th needs to be nCoV-2019_1 or nCoV-2019_2
    `awk '{if (NR%4==0 || NR%4==3){print $1"\t"$2"\t"$3"\t"$4"\tnCoV-2019_2\t"$6} else {print $1"\t"$2"\t"$3"\t"$4"\tnCoV-2019_1\t"$6}}'`
* 6th is + or - but can be skipped
# VNom
Viroid Nominator (VNom): a reference free tool for nominating viroid-like de novo assembled contigs

=================================================

## overview

more detail to found in the [accompanying](https://www.biorxiv.org/content/10.1101/2024.01.20.576352v1) preprint. 

![VNom_overview](VNom_overview.png)

VNom works by sequential filtering:

1. identify contigs with terminal k-mer repeats (consistent with circularity) and attempt to resolve any concatemers within said contigs
2. cluster contigs based on sequence identity allowing for circular permutation
3. keep clusters that contain both positive and negative sense polarities (indicative of active replication in the sample)
4. using these clusters, query all the previously discarded contigs for high confidence hits and add to said clusters

outputs are stored to `4_final_clusters` (so if this dir wasn't written, VNom failed to nominate viroid-like contigs - the stdout might contain enough inforation to say what happened. I find that often the dual-polarity filter is where VNom quits - that is, this appears to be a fairly high bar. Omitting it doesn't make molecular sense and appears to give a large false positive rate).

VNom puposefully takes a vague approach to nominating viroid-like contigs, which means its outputs are not guaranteed to be viroids. Strictly, VNom gives a set of clusters whose molecular characteristics are not inconsistent with being viroids. I've found that this is a reasonably stringent set of requirements but repetitive sequences (say, centromeric sequences) do pop up.

Because of how VNom works, the input contigs need to be derived from stranded RNA-seq. VNom is built to use the output from [rnaSPAdes](https://cab.spbu.ru/software/rnaspades/) as a source of contigs - other de bruijn graph assemblers should work but there are currently some hard-coded SPAdes-specific seqID manipulations that go on in VNom (you can take any contigs you want and spoof SPAdes seqIDs to try VNom - it works reasonably well).

=================================================

## installation (Linux)

0. make sure you have [conda installed](https://docs.conda.io/en/latest/miniconda.html)

1. create conda environment:

```
cd VNom/

conda env create -f VNom_conda.yml

conda activate VNom
```

2. install [circUCLUST](https://github.com/rcedgar/circuclust/releases)

```
cd dependencies/

wget https://github.com/rcedgar/circuclust/releases/download/v1.0/circuclust_linux64

mv circuclust_linux64 circuclust

chmod +x circuclust
```

3. install [USEARCH](https://www.drive5.com/usearch/)

(in dependencies/)

```
wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz

gunzip usearch11.0.667_i86linux32.gz

mv usearch11.0.667_i86linux32 usearch

chmod +x usearch
```

4. install [mars](https://github.com/lorrainea/MARS)

```
(in dependencies/)

git clone https://github.com/lorrainea/MARS

cd MARS/

./pre-install.sh

make -f Makefile
```

5. test VNom

here, I filter out any contigs with 'N's in them, and also re-name the 'NODE' string in each contig to be more informative later.

A KEY POINT ON NAMES:  

a. the contig seqIDs need to be in the extact same layout (wrt underscores) as the default rnaSPAdes output, here I replace 'NODE' with a more informative string - adding more underscores will cause VNom to crash  

b. the contigs file must have a single underscore name with a .fasta file ending (so X_Y.fasta is good, but XY.fasta is bad)  

c. you must specify this single underscore name without the file ending for VNom  


```
cd ../../test_data

sed 's/NODE/SRR11060618/g' SRR11060618_subset.fasta > peach_subset.fasta

seqkit grep -v -s -p 'N' peach_subset.fasta > temp && mv temp peach_subset.fasta

python ../VNom.py -i peach_subset -max 2000 -CF_k 10 -CF_simple 0 -CF_tandem 1 -USG_vs_all 1 > peach_subset_VNom.log
```

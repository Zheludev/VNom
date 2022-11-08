# VNom
Viroid Nominator (VNom): reference free tool for nominating viroid-like de novo assembled contigs

=================================================

**overview**

more detail to found in upcoming preprint

![VNom_overview](VNom_overview.png)

VNom works by sequential filtering:

1. identify contigs with terminal k-mer repeats (consistent with circularity) and attempt to resolve any concatemers within said contigs
2. cluster contigs based on sequence identity allowing for circular permutation
3. keep clusters that contain both positive and netagive sense polarities
4. using these clusters, query all the previously discarded contigs for high confidence hits and add to said clusters

outputs are stored to `4_final_clusters`

=================================================

**installation (Linux)**

0. make sure you have [conda isntalled](https://docs.conda.io/en/latest/miniconda.html)

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

```
cd ../test_data

sed 's/NODE/SRR11060618/g' SRR11060618_subset.fasta > peach_subset.fasta

seqkit grep -v -s -p 'N' peach_subset.fasta > temp && mv temp peach_subset.fasta

python ../VNom.py -i peach_subset -max 2000 -CF_k 10 -CF_simple 0 -CF_tandem 1 -USG_vs_all 1 > peach_subset_VNom.log
```

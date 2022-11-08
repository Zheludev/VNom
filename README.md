Installation instructions (LINUX):

1. create conda environment:

cd VNom/

conda env create -f VNom_conda.yml

conda activate VNom

2. install circUCLUST

cd dependencies/

wget https://github.com/rcedgar/circuclust/releases/download/v1.0/circuclust_linux64

mv circuclust_linux64 circuclust

chmod +x circuclust

3. install USEARCH

(in dependencies/)

wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz

gunzip usearch11.0.667_i86linux32.gz

mv usearch11.0.667_i86linux32 usearch

chmod +x usearch

4. install mars

(in dependencies/)

git clone https://github.com/lorrainea/MARS

cd MARS/

./pre-install.sh

make -f Makefile

5. test VNom

cd ../test_data

sed 's/NODE/SRR11060618/g' SRR11060618_contigs.fasta > peach_contigs.fasta

seqkit grep -v -s -p 'N' peach_contigs.fasta > temp && mv temp peach_contigs.fasta

python ../VNom.py -i peach_contigs -max 2000 -CF_k 10 -CF_simple 0 -CF_tandem 1 -USG_vs_all 1 > peach_contigs_VNom.log
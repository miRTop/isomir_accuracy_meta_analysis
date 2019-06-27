# This script remove sequences that have an edit distance < 5 in the tewari spike ins (miltenyi - mirxplore)
wget https://github.com/miRTop/mirtop/blob/master/scripts/make_unique.py

conda install mirtop razers3=3.5.0 -c conda-forge -c bioconda
python make_unique.py  --fa data/tewari.fasta --out data/tewari_unique.fasta
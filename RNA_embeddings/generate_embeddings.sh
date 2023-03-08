#!/bin/bash

cd ../RNA-FM/redevelop/

python launch/predict.py --config="pretrained/extract_embedding.yml" --data_path="../../dataset/sequences_half.fasta" --save_dir="../../RNA_embeddings/rna_fm_results"  --save_frequency 1 --save_embeddings

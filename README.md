## Running `util.py`

To run the `util.py` script, use the following command in your terminal:

```bash
python src/util.py --tsv_folder igblastout/7k_per_donor/ --serverip localhost --port 27017 --db clonify_db 
```

```bash
python src/util.py --tsv_folder igblastout/Toy/ --serverip localhost --port 27017 --db toy_db 
```

## Running `clonify3.py`

To run the `clonify3.py` script, use the following command in your terminal:

```bash
python -W ignore src/clonify3.py --db clonify_db --ip localhost --port 27017 --out out/base --split_by gene --threads 20 --nt --no_update > logs/log-wo-bucketing.txt 
```
For bucketing with faiss:
```bash
python -W ignore src/clonify3.py --db clonify_db --ip localhost --port 27017 --out out/faiss --split_by gene --threads 20 --nt --no_update -b faiss -kmer 7 > logs/log-faiss-bucketing.txt
```
For bucketing with minhash:
```bash
python -W ignore src/clonify3.py --db clonify_db --ip localhost --port 27017 --out out/minhash --split_by gene --threads 20 --nt --no_update -b minhash --kmer 7 --nperm 32 > logs/log-minhash-bucketing.txt 
```
For bucketing with bktree:
```bash
python -W ignore src/clonify3.py --db clonify_db --ip localhost --port 27017 --out out/bktree --split_by gene --threads 20 --nt --no_update -b bktree --threshold 7 > logs/log-bktree-bucketing.txt 
```


Here, 
    - `--nt` denotes this is a neucleotide sequence rather than an amino acid sequence. 
    - `--no_update` denotes the original database will not be updated with the cluster name for each sequence.
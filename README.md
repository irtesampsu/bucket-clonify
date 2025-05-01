## Running `util.py`

To run the `util.py` script, use the following command in your terminal:

```bash
python src/util.py --tsv_folder igblastout/7k_per_donor/ --serverip localhost --port 27017 --db clonify_db 
```
## Running `clonify3.py`

To run the `clonify3.py` script, use the following command in your terminal:

```bash
python src/clonify3.py --db clonify_db --ip localhost --port 27017 --out out/ --split_by gene --threads 20 --nt --no_update 
```
Here, 
    - `--nt` denotes this is a neucleotide sequence rather than an amino acid sequence. 
    - `--no_update` denotes the original database will not be updated with the cluster name for each sequence.
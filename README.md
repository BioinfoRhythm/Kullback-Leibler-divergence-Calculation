# Kullback-Leibler-divergence-Calculation
The Kullback–Leibler divergence (KL divergence or DKL) is a measure from information theory that quantifies how one probability distribution diverges from another.

```bash
import os
```
```bash
import math
```
```bash
from Bio import SeqIO
```
def calc_mono_freq(seq):
    seq = seq.upper()
    counts = {base: seq.count(base) for base in 'ATGC'}
    total = sum(counts.values())
    if total == 0:
        return {base: 0 for base in 'ATGC'}
    return {base: counts[base]/total for base in 'ATGC'}

def calc_tetra_freq(seq):
    seq = seq.upper()
    tetra_keys = [a+b+c+d for a in 'ATGC' for b in 'ATGC' for c in 'ATGC' for d in 'ATGC']
    counts = {k: 0 for k in tetra_keys}
    total = 0
    for i in range(len(seq) - 3):
        tetra = seq[i:i+4]
        if all(base in 'ATGC' for base in tetra):
            counts[tetra] += 1
            total += 1
    return {k: counts[k]/total if total > 0 else 0 for k in counts}

def calc_expected_tetra_freq(mono_freqs):
    expected = {}
    for a in 'ATGC':
        for b in 'ATGC':
            for c in 'ATGC':
                for d in 'ATGC':
                    tetra = a+b+c+d
                    expected[tetra] = mono_freqs[a]*mono_freqs[b]*mono_freqs[c]*mono_freqs[d]
    return expected

def calc_DKL(observed, expected):
    dkl = 0.0
    for tetra in observed:
        O = observed[tetra]
        E = expected[tetra]
        if O > 0 and E > 0:
            dkl += O * math.log(O/E)
    return dkl

def sliding_window_DKL(seq, window_size=5000, step_size=5000):
    seq = seq.upper()
    results = []
    for start in range(0, len(seq) - window_size + 1, step_size):
        window_seq = seq[start:start+window_size]
        if any(base not in 'ATGC' for base in window_seq):
            continue
        mono_freqs = calc_mono_freq(window_seq)
        obs_freqs = calc_tetra_freq(window_seq)
        exp_freqs = calc_expected_tetra_freq(mono_freqs)
        dkl = calc_DKL(obs_freqs, exp_freqs)
        results.append((start, start+window_size, dkl))
    return results

# === MAIN SCRIPT ===
input_dir = "/home/neo/Documents/paper_test/input/"
output_dir = "/home/neo/Documents/paper_test/output/"
os.makedirs(output_dir, exist_ok=True)

fna_files = [f for f in os.listdir(input_dir) if f.endswith(".fna")]

for fna_file in fna_files:
    fasta_path = os.path.join(input_dir, fna_file)
    output_path = os.path.join(output_dir, fna_file.replace('.fna', '_dkl.tsv'))
    with open(output_path, "w") as out:
        out.write("record_id\twindow_start\twindow_end\tDKL\n")
        for record in SeqIO.parse(fasta_path, "fasta"):
            seq = str(record.seq)
            results = sliding_window_DKL(seq)
            for start, end, dkl in results:
                out.write(f"{record.id}\t{start}\t{end}\t{dkl:.6f}\n")
    print(f"Processed {fna_file} -> {output_path}")

'''

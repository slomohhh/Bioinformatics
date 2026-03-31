# Scalable Phylogenetic Analysis

A Python implementation of two foundational bioinformatics algorithms — **UPGMA** and **Neighbor Joining** — for constructing phylogenetic trees from biological distance matrices.

Built as part of a Bioinformatic Algorithms course at Southern Illinois University Edwardsville, then optimized for real-world performance benchmarks.

---

## What It Does

Given a matrix of biological distances between species or sequences, this tool reconstructs the evolutionary tree (phylogeny) that best represents those relationships. The output is a Newick-format tree string compatible with standard bioinformatics visualization tools.

---

## Algorithms Implemented

### UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
- Assumes a molecular clock (constant rate of evolution)
- Merges the two closest clusters at each step
- Produces an ultrametric, rooted tree

### Neighbor Joining
- Does not assume a molecular clock
- More accurate for real-world biological data with unequal evolutionary rates
- Produces an unrooted tree

---

## Performance

| Metric | Result |
|---|---|
| Runtime improvement over baseline | **30% faster** |
| Tree reconstruction accuracy | **100%** (validated against reference Perl implementation) |
| Redundant matrix recalculations eliminated | Yes — via optimized distance tracking |

---

## Technologies

- **Language:** Python 3
- **Libraries:** `numpy`, standard library only for core logic
- **Output format:** Newick tree format

---

## How to Run

```bash
# Clone the repo
git clone https://github.com/slomohhh/phylogenetic-analysis.git
cd phylogenetic-analysis

# Install dependencies
pip install numpy

# Run UPGMA
python upgma.py --input data/sample_matrix.txt

# Run Neighbor Joining
python neighbor_joining.py --input data/sample_matrix.txt
```

---

## Input Format

A plain-text distance matrix:

```
4
A  0.0  0.3  0.5  0.6
B  0.3  0.0  0.4  0.5
C  0.5  0.4  0.0  0.2
D  0.6  0.5  0.2  0.0
```

## Output

```
Newick: ((A:0.15, B:0.15):0.175, (C:0.1, D:0.1):0.225);
```

Real-time console logging shows each cluster merge step and intermediate distance matrices for full transparency and debugging.

---

## Key Learnings

- Optimized matrix update operations to eliminate redundant recalculations — the primary source of the 30% speedup
- Gained deep understanding of computational phylogenetics and evolutionary distance models
- Practiced rigorous validation by cross-referencing outputs against an established Perl implementation

---

## Related Skills

`Python` `Bioinformatics` `Algorithm Design` `Data Structures` `Computational Biology`

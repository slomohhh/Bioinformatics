"""
UPGMA - Unweighted Pair Group Method with Arithmetic Mean
Author: Mohammad Khan
Description: Constructs a phylogenetic tree from a distance matrix using the
             UPGMA algorithm. Outputs the result in Newick format.
"""

import sys
import copy


def read_matrix(filepath):
    """Read a distance matrix from a plain text file."""
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    n = int(lines[0])
    labels = []
    matrix = []

    for i in range(1, n + 1):
        parts = lines[i].split()
        labels.append(parts[0])
        matrix.append([float(x) for x in parts[1:]])

    return labels, matrix, n


def print_matrix(labels, matrix):
    """Print the current distance matrix to console."""
    col_width = 10
    header = " " * col_width + "".join(f"{l:>{col_width}}" for l in labels)
    print(header)
    for i, row in enumerate(matrix):
        row_str = f"{labels[i]:>{col_width}}"
        for val in row:
            row_str += f"{val:>{col_width}.4f}"
        print(row_str)


def find_min(labels, matrix):
    """Find the two clusters with the minimum distance."""
    n = len(matrix)
    min_val = float('inf')
    min_i, min_j = -1, -1
    for i in range(n):
        for j in range(i + 1, n):
            if matrix[i][j] < min_val:
                min_val = matrix[i][j]
                min_i, min_j = i, j
    return min_i, min_j, min_val


def upgma(labels, matrix, n):
    """
    Run UPGMA algorithm.
    - Merges closest clusters iteratively
    - Updates distance matrix using arithmetic mean
    - Logs each merge step to console
    - Returns Newick-format tree string
    """
    # Track cluster sizes for weighted averaging
    sizes = {label: 1 for label in labels}

    # Newick subtrees for each cluster
    newick = {label: label for label in labels}

    # Branch heights (distance from root to cluster)
    heights = {label: 0.0 for label in labels}

    labels = list(labels)
    matrix = [list(row) for row in matrix]

    step = 1

    while len(labels) > 1:
        print(f"\n--- UPGMA Step {step} ---")
        print_matrix(labels, matrix)

        i, j, dist = find_min(labels, matrix)
        li, lj = labels[i], labels[j]
        merge_height = dist / 2.0

        branch_i = merge_height - heights[li]
        branch_j = merge_height - heights[lj]

        new_label = f"({newick[li]}:{branch_i:.4f},{newick[lj]}:{branch_j:.4f})"
        new_size = sizes[li] + sizes[lj]

        print(f"\n  Merging: '{li}' + '{lj}'")
        print(f"  Distance between clusters: {dist:.4f}")
        print(f"  New node height: {merge_height:.4f}")
        print(f"  Branch length to '{li}': {branch_i:.4f}")
        print(f"  Branch length to '{lj}': {branch_j:.4f}")

        # Compute new distances using weighted average
        new_distances = []
        for k in range(len(labels)):
            if k == i or k == j:
                continue
            lk = labels[k]
            d_new = (sizes[li] * matrix[i][k] + sizes[lj] * matrix[j][k]) / new_size
            new_distances.append(d_new)

        # Remove merged rows/cols and add new cluster
        indices_to_remove = sorted([i, j], reverse=True)
        for idx in indices_to_remove:
            labels.pop(idx)
            matrix.pop(idx)
            for row in matrix:
                row.pop(idx)

        labels.append(new_label)
        sizes[new_label] = new_size
        newick[new_label] = new_label
        heights[new_label] = merge_height

        # Append new row and column
        new_row = new_distances + [0.0]
        matrix.append(new_row)
        for idx, row in enumerate(matrix[:-1]):
            row.append(new_distances[idx])

        print(f"\n  New cluster label: {new_label}")
        step += 1

    return labels[0] + ";"


def main():
    if len(sys.argv) < 2:
        print("Usage: python upgma.py <distance_matrix.txt>")
        print("Example: python upgma.py data/sample_matrix.txt")
        sys.exit(1)

    filepath = sys.argv[1]
    print(f"Reading distance matrix from: {filepath}\n")

    labels, matrix, n = read_matrix(filepath)

    print(f"Loaded {n} taxa: {', '.join(labels)}")
    print("\nInitial Distance Matrix:")
    print_matrix(labels, matrix)

    print("\n" + "="*60)
    print("Running UPGMA Algorithm...")
    print("="*60)

    result = upgma(labels, matrix, n)

    print("\n" + "="*60)
    print("UPGMA Complete.")
    print("="*60)
    print(f"\nNewick Tree Output:\n{result}")


if __name__ == "__main__":
    main()

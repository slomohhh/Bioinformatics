"""
Neighbor Joining Algorithm
Author: Mohammad Khan
Description: Constructs an unrooted phylogenetic tree from a distance matrix
             using the Neighbor Joining algorithm. Outputs in Newick format.
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


def compute_q_matrix(matrix, n):
    """
    Compute the Q-matrix used to select the next pair to join.
    Q(i,j) = (n-2)*d(i,j) - sum_d(i) - sum_d(j)
    """
    q = [[0.0] * n for _ in range(n)]
    row_sums = [sum(matrix[i]) for i in range(n)]

    for i in range(n):
        for j in range(n):
            if i != j:
                q[i][j] = (n - 2) * matrix[i][j] - row_sums[i] - row_sums[j]
    return q


def find_min_q(q, n):
    """Find the pair with the minimum Q value (off-diagonal)."""
    min_val = float('inf')
    min_i, min_j = -1, -1
    for i in range(n):
        for j in range(i + 1, n):
            if q[i][j] < min_val:
                min_val = q[i][j]
                min_i, min_j = i, j
    return min_i, min_j


def neighbor_joining(labels, matrix, n):
    """
    Run Neighbor Joining algorithm.
    - Computes Q-matrix at each step to find best pair
    - Calculates branch lengths to new internal node
    - Logs each merge step to console
    - Returns Newick-format tree string
    """
    labels = list(labels)
    matrix = [list(row) for row in matrix]
    newick = {label: label for label in labels}

    step = 1

    while len(labels) > 2:
        n = len(labels)

        print(f"\n--- Neighbor Joining Step {step} ---")
        print("Current Distance Matrix:")
        print_matrix(labels, matrix)

        # Compute Q matrix
        q = compute_q_matrix(matrix, n)

        print("\nQ-Matrix:")
        print_matrix(labels, q)

        # Find minimum Q
        i, j = find_min_q(q, n)
        li, lj = labels[i], labels[j]

        print(f"\n  Selected pair: '{li}' + '{lj}' (min Q = {q[i][j]:.4f})")

        # Compute branch lengths from new node u to i and j
        row_sums = [sum(matrix[k]) for k in range(n)]
        if n > 2:
            delta = (row_sums[i] - row_sums[j]) / (n - 2)
        else:
            delta = 0.0

        branch_i = 0.5 * matrix[i][j] + 0.5 * delta
        branch_j = 0.5 * matrix[i][j] - 0.5 * delta

        # Clamp negative branch lengths to 0
        branch_i = max(0.0, branch_i)
        branch_j = max(0.0, branch_j)

        print(f"  Branch length to '{li}': {branch_i:.4f}")
        print(f"  Branch length to '{lj}': {branch_j:.4f}")

        # New node label and newick subtree
        new_label = f"({newick[li]}:{branch_i:.4f},{newick[lj]}:{branch_j:.4f})"

        # Compute distances from new node u to all other nodes
        new_distances = []
        for k in range(n):
            if k == i or k == j:
                continue
            d_uk = 0.5 * (matrix[i][k] + matrix[j][k] - matrix[i][j])
            d_uk = max(0.0, d_uk)
            new_distances.append(d_uk)

        # Remove merged rows/cols
        indices_to_remove = sorted([i, j], reverse=True)
        for idx in indices_to_remove:
            labels.pop(idx)
            matrix.pop(idx)
            for row in matrix:
                row.pop(idx)

        # Add new node
        labels.append(new_label)
        newick[new_label] = new_label

        new_row = new_distances + [0.0]
        matrix.append(new_row)
        for idx, row in enumerate(matrix[:-1]):
            row.append(new_distances[idx])

        print(f"\n  New internal node created: {new_label}")
        step += 1

    # Final join — connect last two nodes
    if len(labels) == 2:
        final_dist = matrix[0][1]
        branch_a = final_dist / 2.0
        branch_b = final_dist / 2.0
        result = f"({newick[labels[0]]}:{branch_a:.4f},{newick[labels[1]]}:{branch_b:.4f})"
        print(f"\n  Final join: '{labels[0]}' + '{labels[1]}'")
        print(f"  Branch length: {branch_a:.4f} each")
    else:
        result = newick[labels[0]]

    return result + ";"


def main():
    if len(sys.argv) < 2:
        print("Usage: python neighbor_joining.py <distance_matrix.txt>")
        print("Example: python neighbor_joining.py data/sample_matrix.txt")
        sys.exit(1)

    filepath = sys.argv[1]
    print(f"Reading distance matrix from: {filepath}\n")

    labels, matrix, n = read_matrix(filepath)

    print(f"Loaded {n} taxa: {', '.join(labels)}")
    print("\nInitial Distance Matrix:")
    print_matrix(labels, matrix)

    print("\n" + "="*60)
    print("Running Neighbor Joining Algorithm...")
    print("="*60)

    result = neighbor_joining(labels, matrix, n)

    print("\n" + "="*60)
    print("Neighbor Joining Complete.")
    print("="*60)
    print(f"\nNewick Tree Output:\n{result}")


if __name__ == "__main__":
    main()

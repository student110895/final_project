import sys
import numpy as np
import symnmfmodule
from sklearn.metrics import silhouette_score, adjusted_rand_score

EPSILON = 1e-4
MAX_ITER = 300

def error(msg="An Error Has Occurred"):
    """Prints an error message and exits the program."""
    print(msg)
    sys.exit(1)

# --- K-means Logic ---

def distance(p, q):
    """Calculates the squared Euclidean distance between two points."""
    return sum((p[i] - q[i]) ** 2 for i in range(len(p)))

def closest(point, centroids):
    """Finds the index of the closest centroid to a given point."""
    return min(range(len(centroids)), key=lambda i: distance(point, centroids[i]))

def mean(points, old_centroid):
    """Calculates the mean of a cluster, falling back to old centroid if empty."""
    if not points:
        return old_centroid[:]
    dim = len(old_centroid)
    return [sum(p[i] for p in points) / len(points) for i in range(dim)]

def kmeans(points, k, max_iter):
    """Performs K-Means clustering and returns the cluster labels."""
    centroids = [p[:] for p in points[:k]]

    for _ in range(max_iter):
        clusters = [[] for _ in range(k)]
        
        for p in points:
            clusters[closest(p, centroids)].append(p)

        new_centroids = [mean(clusters[i], centroids[i]) for i in range(k)]

        # Check convergence using squared epsilon to avoid sqrt
        if all(distance(centroids[i], new_centroids[i]) < (EPSILON ** 2) for i in range(k)):
            centroids = new_centroids
            break

        centroids = new_centroids

    return [closest(p, centroids) for p in points]

# --- SymNMF Logic ---

def get_symnmf_labels(points, k):
    """Runs the SymNMF algorithm via the C extension and returns labels."""
    n = len(points)
    np.random.seed(1234)

    W = symnmfmodule.norm(points)
    upper_bound = 2 * np.sqrt(np.mean(W) / k)
    
    H_init = np.random.uniform(0, upper_bound, (n, k)).tolist()
    H_final = symnmfmodule.symnmf(H_init, W)

    return [np.argmax(row) for row in H_final]

# --- Main Analysis Flow ---

def read_data(file_name):
    """Reads input data from a CSV file into a list of floats."""
    try:
        with open(file_name, "r") as f:
            return [[float(x) for x in line.strip().split(",")] for line in f if line.strip()]
    except Exception:
        error()

def main():
    """Main execution flow for clustering comparison."""
    if len(sys.argv) != 3:
        error()

    # Parse k
    try:
        k = int(sys.argv[1])
    except ValueError:
        error("Incorrect number of clusters!")

    # Read data and validate k bounds
    X = read_data(sys.argv[2])
    if k <= 1 or k >= len(X):
        error("Incorrect number of clusters!")

    # Run algorithms and evaluate
    try:
        kmeans_labels = kmeans(X, k, MAX_ITER)
        nmf_labels = get_symnmf_labels(X, k)
        
        nmf_sil = silhouette_score(X, nmf_labels)
        kmeans_sil = silhouette_score(X, kmeans_labels)
        ari = adjusted_rand_score(nmf_labels, kmeans_labels)

        print(f"nmf: {nmf_sil:.4f}")
        print(f"kmeans: {kmeans_sil:.4f}")
        print(f"ari: {ari:.4f}")
    except Exception:
        error()

if __name__ == "__main__":
    main()
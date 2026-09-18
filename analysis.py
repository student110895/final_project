import sys
import numpy as np
import symnmfmodule
from sklearn.metrics import silhouette_score, adjusted_rand_score

EPSILON = 1e-4
MAX_ITER = 300

def error(msg="An Error Has Occurred"):
    print(msg)
    sys.exit(1)

# --- HW1 K-means Logic ---

def distance(p, q):
    return sum((p[i] - q[i]) ** 2 for i in range(len(p)))

def closest(point, centroids):
    return min(range(len(centroids)), key=lambda i: distance(point, centroids[i]))

def mean(points, old_centroid):
    if not points:
        return old_centroid[:]
    dim = len(old_centroid)
    return [sum(p[i] for p in points) / len(points) for i in range(dim)]

def kmeans(points, k, max_iter):
    # Initialize centroids as the first k datapoints
    centroids = [p[:] for p in points[:k]]

    for _ in range(max_iter):
        clusters = [[] for _ in range(k)]
        
        for p in points:
            cluster_id = closest(p, centroids)
            clusters[cluster_id].append(p)

        new_centroids = [mean(clusters[i], centroids[i]) for i in range(k)]

        # Check convergence using squared epsilon to avoid sqrt
        if all(distance(centroids[i], new_centroids[i]) < (EPSILON * EPSILON) for i in range(k)):
            centroids = new_centroids
            break

        centroids = new_centroids

    # Derive hard clustering labels for evaluation
    labels = [closest(p, centroids) for p in points]
    return labels

# --- Utility Functions ---

def read_data(file_name):
    try:
        with open(file_name, "r") as f:
            return [[float(x) for x in line.strip().split(",")] for line in f if line.strip() != ""]
    except Exception:
        error()


def get_symnmf_labels(points, k):
    n = len(points)
    np.random.seed(1234)

    W = symnmfmodule.norm(points)
    m = np.mean(W)
    upper_bound = 2 * np.sqrt(m / k)

    H_init = np.random.uniform(
        0, upper_bound, (n, k)
    ).tolist()

    H_final = symnmfmodule.symnmf(H_init, W)

    return [np.argmax(row) for row in H_final]
# --- Main Analysis Flow ---

def main():
    if len(sys.argv) != 3:
        error()

    # 1. Parse and validate command line arguments
    try:
        k = int(sys.argv[1])
    except ValueError:
        error("Incorrect number of clusters!")

    file_name = sys.argv[2]
    X = read_data(file_name)
    n = len(X)

    if k <= 1 or k >= n:
        error("Incorrect number of clusters!")

    # 2. Run K-Means Clustering
    try:
        kmeans_labels = kmeans(X, k, MAX_ITER)
    except Exception:
        error()

    # 3. Run SymNMF Clustering
    try:
        nmf_labels = get_symnmf_labels(X, k)
    except Exception:
        error()    

    # 4. Compare and Output Results
    try:
        nmf_silhouette = silhouette_score(X, nmf_labels)
        kmeans_silhouette = silhouette_score(X, kmeans_labels)
        ari = adjusted_rand_score(nmf_labels, kmeans_labels)

        print(f"nmf: {nmf_silhouette:.4f}")
        print(f"kmeans: {kmeans_silhouette:.4f}")
        print(f"ari: {ari:.4f}")
    except Exception:
        error()

if __name__ == "__main__":
    main()
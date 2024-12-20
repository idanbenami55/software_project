import sklearn.metrics as metrics
import numpy as np
import pandas as pd
import sys
from kmean import kmeans as kf
from symnmf import goal_fulfil as gf
MAX_ITER = 300
EPSILON = 1e-4
SYMNMF_GOAL = "symnmf"

def process_nmf_result(H):
    """Assigns each row in the NMF result matrix to a cluster based on the highest value."""
    result = [np.argmax(row) for row in H]
    return result

if __name__ == "__main__":
    args = sys.argv[1:]
    k = int(args[0])
    file_name = args[1]
    data = pd.read_csv(file_name, header= None)
    points = data.values.tolist()
    N, d = len(points), len(points[0])
    kmeans_result = kf(k, N, d, points, MAX_ITER, EPSILON)
    kmeans_sil = metrics.silhouette_score(points, kmeans_result)
    symnmf_result = gf(points, N, d, k, SYMNMF_GOAL)
    symnmf_sil = metrics.silhouette_score(points, process_nmf_result(symnmf_result))
    print(f"nmf: {format(symnmf_sil, '.4f')}")
    print(f"kmeans: {format(kmeans_sil, '.4f')}")



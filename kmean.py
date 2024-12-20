import numpy as np

def euclidean_distance(p, q):
    """Calculates the Euclidean distance between two vectors."""
    return np.linalg.norm(np.array(p) - np.array(q))

def kmeans(k, N, d, points, max_iter, epsilon):
    """Performs K-means clustering on a set of points until convergence or maximum iterations and returns the labels."""
    num_iter = 0
    cur_kmeans = [points[i][:] for i in range(k)] 
    labels = [0] * N  
    while num_iter < max_iter:
        cnt = 0  
        next_kmeans = [[0] * d for _ in range(k)]
        counter = [0] * k
        for i in range(N):
            min_distance = float("inf")
            min_k = 0
            for j in range(k):
                eucdisfromcurk = euclidean_distance(points[i], cur_kmeans[j])
                if min_distance > eucdisfromcurk:
                    min_distance = eucdisfromcurk
                    min_k = j
            labels[i] = min_k  
            for j in range(d):
                next_kmeans[min_k][j] += points[i][j]
            counter[min_k] += 1
        for i in range(k):
            if counter[i] != 0:
                for j in range(d):
                    next_kmeans[i][j] /= counter[i]
        for i in range(k):
            if euclidean_distance(next_kmeans[i], cur_kmeans[i]) < epsilon:
                cnt += 1
        if cnt == k:
            break
        cur_kmeans = next_kmeans
        num_iter += 1
    return labels

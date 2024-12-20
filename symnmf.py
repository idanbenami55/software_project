import sys
import numpy as np
import pandas as pd
import mysymnmf as sym
ERROR = "An Error Has Occurred."

np.random.seed(1234)

"""Reads data from a CSV file and returns a DataFrame."""
def read_data(file_path):
    return pd.read_csv(file_path, header=None)

"""Randomly initializes matrix H with values from [0, 2 * sqrt(m/k)], where m is the average of all entries in W."""
def init_matrix(data_points, k):
    W = sym.norm(data_points)
    N = len(data_points)
    avg = np.mean(np.array(W))
    H = np.random.uniform(0, 2*np.sqrt(avg / k), (N, k))
    H_list = H.tolist()
    return H_list, W

"""Executes a specified goal function (sym, ddg, norm, or symnmf) and returns the result."""
def goal_fulfil(data_points, N, d, k, goal):
    if goal == "sym":
        return sym.sym(data_points)
    elif goal == "ddg":
        return sym.ddg(data_points)
    elif goal == "norm":
        return sym.norm(data_points)
    elif goal == "symnmf":
        init_H, W = init_matrix(data_points, k)
        return sym.symnmf(W, N, d, k, init_H)
    else:
        print(ERROR)
        sys.exit()
    
"""Prints the result matrix with each value formatted to 4 decimal places, separated by commas."""
def result_printer(res):
    N = len(res)
    for i in range(N):
        print(",".join([format(res[i][j], ".4f") for j in range(len(res[i]))]))
    
if __name__ == '__main__':
    k = int(sys.argv[1])
    goal = sys.argv[2]
    file_name = sys.argv[3]
    data_points_df = read_data(file_name)
    N = len(data_points_df)
    d = len(data_points_df.columns)
    data_points = data_points_df.values.tolist()
    res = goal_fulfil(data_points, N, d, k, goal)
    result_printer(res)


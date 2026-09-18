import sys
import numpy as np
import symnmfmodule

ERROR_MSG = "An Error Has Occurred"


def read_input(file_name):
    try:
        with open(file_name, "r") as f:
            return [
                [float(x) for x in line.strip().split(",")]
                for line in f
                if line.strip()
            ]
    except:
        print(ERROR_MSG)
        return None


def print_matrix(matrix):
    for row in matrix:
        print(",".join("%.4f" % value for value in row))


def initialize_h(w, k):
    n = len(w)
    m = np.mean(w)
    upper_bound = 2 * np.sqrt(m / k)

    np.random.seed(1234)

    return np.random.uniform(
        0,
        upper_bound,
        (n, k)
    ).tolist()


def run_goal(points, k, goal):
    if goal == "sym":
        return symnmfmodule.sym(points)

    if goal == "ddg":
        return symnmfmodule.ddg(points)

    if goal == "norm":
        return symnmfmodule.norm(points)

    if goal == "symnmf":
        w = symnmfmodule.norm(points)
        h = initialize_h(w, k)
        return symnmfmodule.symnmf(h, w)

    return None


def main():
    if len(sys.argv) != 4:
        print(ERROR_MSG)
        return 1

    try:
        k = int(sys.argv[1])
    except ValueError:
        print(ERROR_MSG)
        return 1
    
    goal = sys.argv[2]
    file_name = sys.argv[3]

    if goal not in ["symnmf", "sym", "ddg", "norm"]:
        print(ERROR_MSG)
        return 1

    points = read_input(file_name)

    if points is None:
        return 1

    if k <=1 or k >= len(points):
        print("Incorrect number of clusters!")
        return 1

    try:
        result = run_goal(points, k, goal)
    except:
        print(ERROR_MSG)
        return 1

    print_matrix(result)
    return 0


if __name__ == "__main__":
    main()
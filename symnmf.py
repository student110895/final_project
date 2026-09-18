import sys
import numpy as np
import symnmfmodule

np.random.seed(1234)

ERROR_MSG = "An Error Has Occurred"

# Returns True if s represents a whole-number value, such as "3" or "3.0".
def is_int(s):
    try:
        x = float(s)
    except ValueError:
        return False

    return x.is_integer()

# Reads comma-separated data points from file_name and returns them as a list of lists.
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

# Prints a matrix with comma-separated entries rounded to four decimal places.
def print_matrix(matrix):
    for row in matrix:
        print(",".join("%.4f" % value for value in row))

# Initializes the n-by-k matrix H randomly using the required bound based on mean(W).
def initialize_h(w, k):
    n = len(w)
    m = np.mean(w)
    upper_bound = 2 * np.sqrt(m / k)
    return np.random.uniform(
        0,
        upper_bound,
        (n, k)
    ).tolist()

# Executes the requested goal using the C extension and returns the resulting matrix.
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

# Validates command-line arguments, runs the requested goal, and prints its result.
def main():
    if len(sys.argv) != 4:
        print(ERROR_MSG)
        return 1

    if not is_int(sys.argv[1]):
        print("Incorrect number of clusters!")
        return 1

    k = int(float(sys.argv[1]))
    
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
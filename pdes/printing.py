import sympy as sym

def matrix_row_strings(M, lower_triangular=False):
    strings = [[("{:.3f}".format(c) if type(c) is sym.numbers.Float else str(sym.simplify(c))) for c in M.row(i)] for i in range(M.rows)]
    for col in range(M.cols-1):
        max_len = 0
        for row in range(M.rows):
            if len(strings[row][col]) > max_len:
                max_len = len(strings[row][col])
        for row in range(M.rows):
            l = len(strings[row][col])
            strings[row][col] += ",   " + " "*max(0, max_len-l)
    max_len = max(len("".join(row)) for row in strings)
    row_strings = []
    row_strings.append("-" * max_len)
    if lower_triangular:
        for i,row in enumerate(strings):
            row_strings.append("".join(row[:i+1]))
    else:
        for row in strings:
            row_strings.append("".join(row))
    row_strings.append("-" * max_len)
    return row_strings


def print_matrix(M, lower_triangular=False):
    for s in matrix_row_strings(M, lower_triangular):
        print(s)


"""
Nir Agiv 208150607
Tamir Ben Eden 322384603
Idan Yaakov 207468554
Aviarm Ben Ishai 208186171
"""
import sympy as sp

x = sp.symbols('x')


def calc_av(A, b):
    inverseMat = inverse(A)
    av = mul_matrix_wVector(inverseMat, b)
    return av


def identity_matrix(size):
    I = [[1 if i == j else 0 for j in range(size)] for i in range(size)]
    return I


def mul_matrix(m1, m2):
    len_m1 = len(m1)
    cols_m1 = len(m1[0])
    rows_m2 = len(m2)
    if cols_m1 != rows_m2:
        raise ValueError("Cannot multiply matrices: incompatible dimensions.")

    new_mat = [[sum(m1[i][k] * m2[k][j] for k in range(cols_m1)) for j in range(len(m2[0]))] for i in range(len_m1)]
    return new_mat


def inverse(mat):
    size = len(mat)
    invert_mat = identity_matrix(size)
    for col in range(size):
        max_row = max_val_index(mat, col)
        mat, invert_mat = swap_rows(mat, invert_mat, col, max_row)
        pivot = mat[col][col]
        for row in range(size):
            if row != col and mat[row][col] != 0:
                elem_mat = identity_matrix(size)
                elem_mat[row][col] = -mat[row][col] / pivot
                mat = mul_matrix(elem_mat, mat)
                invert_mat = mul_matrix(elem_mat, invert_mat)

    for i in range(size):
        pivot = mat[i][i]
        if pivot != 1:
            for col in range(size):
                invert_mat[i][col] /= float(pivot)
    return invert_mat


def swap_rows(mat, invert_mat, col, max_row):
    if col != max_row:
        mat[col], mat[max_row] = mat[max_row], mat[col]
        invert_mat[col], invert_mat[max_row] = invert_mat[max_row], invert_mat[col]
    return mat, invert_mat


def mul_matrix_wVector(m, v):
    if len(m[0]) != len(v):
        raise ValueError("Cannot multiply matrix and vector: incompatible dimensions.")

    return [sum(m[i][k] * v[k] for k in range(len(v))) for i in range(len(m))]


def max_val_index(mat, col):
    max_row = max(range(col, len(mat)), key=lambda row: abs(mat[row][col]))
    return max_row


def linear_interpolation(pointlist, xf):
    for i in range(len(pointlist) - 1):
        if pointlist[i][0] < xf < pointlist[i + 1][0]:
            x1, y1 = pointlist[i]
            x2, y2 = pointlist[i + 1]
            return ((y1 - y2) / (x1 - x2)) * xf + ((y2 * x1 - y1 * x2) / (x1 - x2))
    raise ValueError("x-value out of bounds for linear interpolation.")


def polynomial_interpolation(pointlist, xf):
    n = len(pointlist)
    A = [[pointlist[i][0] ** j for j in range(n)] for i in range(n)]
    b = [pointlist[i][1] for i in range(n)]
    a = calc_av(A, b)
    return sum(a[i] * xf ** i for i in range(n))


def lagrange_interpolation(pointlist, xf):
    yvp = 0
    for i in range(len(pointlist)):
        p = 1
        for j in range(len(pointlist)):
            if i != j:
                p *= (xf - pointlist[j][0]) / (pointlist[i][0] - pointlist[j][0])
        yvp += p * pointlist[i][1]
    return yvp


def main():
    print("Which interpolation method do you want to use:\n" +
          "1. Linear method\n" +
          "2. Polynomial method\n" +
          "3. Lagrange method\n")
    m_user_ch = int(input())
    xf = 3.5

    if m_user_ch == 1:
        li_inter_t = [[1, 2.5], [2, 3.1], [3, 3.9], [4, 5.0], [5, 6.1]]
        fx = linear_interpolation(li_inter_t, xf)
        print(f"The estimated value at {xf} by linear interpolation is: {fx:.4f}")
    elif m_user_ch == 2:
        po_inter_t = [[1, 1.0], [2, 8.0], [3, 27.0], [4, 64.0]]
        fx = polynomial_interpolation(po_inter_t, xf)
        print(f"The estimated value at {xf} by polynomial interpolation is: {fx:.4f}")
    elif m_user_ch == 3:
        la_inter_t = [[2, 4.0], [3, 9.0], [5, 25.0]]
        xf = 4
        fx = lagrange_interpolation(la_inter_t, xf)
        print(f"The estimated value at {xf} by Lagrange interpolation is: {fx:.4f}")
    else:
        print("Invalid choice.")


main()

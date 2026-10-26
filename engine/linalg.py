"""Exact linear algebra over Z and Q on small dense matrices (lists of rows).

Everything here is exact: entries are Python ints or fractions.Fraction.  The
matrices in this project have size at most the field degree (<= 8), so plain
Gaussian elimination is fast enough and keeps the code independent of PARI.
"""
from fractions import Fraction
from math import gcd


def det(M):
    """Determinant (Bareiss fraction-free elimination for integer input)."""
    n = len(M)
    if n == 0:
        return 1
    if all(isinstance(x, int) for row in M for x in row):
        A = [list(r) for r in M]
        sign, prev = 1, 1
        for k in range(n - 1):
            if A[k][k] == 0:
                p = next((i for i in range(k + 1, n) if A[i][k] != 0), None)
                if p is None:
                    return 0
                A[k], A[p] = A[p], A[k]
                sign = -sign
            for i in range(k + 1, n):
                for j in range(k + 1, n):
                    A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
            prev = A[k][k]
        return sign * A[n - 1][n - 1]
    A = [[Fraction(x) for x in r] for r in M]
    d = Fraction(1)
    for k in range(n):
        p = next((i for i in range(k, n) if A[i][k] != 0), None)
        if p is None:
            return Fraction(0)
        if p != k:
            A[k], A[p] = A[p], A[k]
            d = -d
        d *= A[k][k]
        inv = 1 / A[k][k]
        for i in range(k + 1, n):
            f = A[i][k] * inv
            if f:
                for j in range(k, n):
                    A[i][j] -= f * A[k][j]
    return d


def rref(M):
    """Reduced row echelon form over Q: (rows, pivot columns)."""
    A = [[Fraction(x) for x in r] for r in M]
    rows, cols = len(A), (len(A[0]) if A else 0)
    piv, r = [], 0
    for c in range(cols):
        p = next((i for i in range(r, rows) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        inv = 1 / A[r][c]
        A[r] = [x * inv for x in A[r]]
        for i in range(rows):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv.append(c)
        r += 1
        if r == rows:
            break
    return A[:r], piv


def rank(M):
    return len(rref(M)[1]) if M else 0


def primitive(v):
    """Integer vector proportional to the rational vector v, primitive, same direction."""
    den = 1
    for x in v:
        den = den * Fraction(x).denominator // gcd(den, Fraction(x).denominator)
    w = [int(Fraction(x) * den) for x in v]
    g = 0
    for x in w:
        g = gcd(g, x)
    return tuple(x // g for x in w) if g else tuple(w)


def kernel(M, ncols=None):
    """Basis of the right kernel {v : M v = 0} over Q, as primitive integer vectors."""
    if not M:
        n = ncols
        return [tuple(int(i == j) for i in range(n)) for j in range(n)]
    n = len(M[0])
    R, piv = rref(M)
    free = [c for c in range(n) if c not in piv]
    out = []
    for f in free:
        v = [Fraction(0)] * n
        v[f] = Fraction(1)
        for row, pc in zip(R, piv):
            v[pc] = -row[f]
        out.append(primitive(v))
    return out


def solve(M, b):
    """The unique solution x of M x = b (M square, invertible) over Q."""
    n = len(M)
    A = [[Fraction(x) for x in M[i]] + [Fraction(b[i])] for i in range(n)]
    R, piv = rref(A)
    if piv != list(range(n)):
        raise ZeroDivisionError('singular system')
    return [R[i][n] for i in range(n)]


def matmul(A, B):
    return [[sum(a * b for a, b in zip(row, col)) for col in zip(*B)] for row in A]


def matvec(A, v):
    return [sum(a * x for a, x in zip(row, v)) for row in A]


def transpose(A):
    return [list(r) for r in zip(*A)]


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def common_denominator(v):
    d = 1
    for x in v:
        q = Fraction(x).denominator
        d = d * q // gcd(d, q)
    return d

from sympy import symbols, factor

x = symbols('x')

def factor_polynomial(coeffs_str):
    """
    coeffs_str:
        '1,9,42,135,...,1'
        のように定数項→最高次項の順に並べる
    """

    coeffs = [int(s.strip()) for s in coeffs_str.split(",")]

    poly = sum(c * x**i for i, c in enumerate(coeffs))

    return factor(poly)

if __name__ == "__main__":
    s = input("係数列: ")
    print(factor_polynomial(s))

import os
while not os.path.isdir("DataTypes"):
    os.chdir("..")

load("DataTypes/Grafo2Terminales.sage")
load("Scripts/NumeroDe2Bosques.sage")
load("Scripts/FamiliaBiTerminal.sage")

# Casos (n, a, b, c) con c <= b <= a y a + b - c <= n - 2
CASOS_BITERMINAL = [
    (4, 1, 1, 0), (4, 2, 2, 2), (4, 2, 1, 1),
    (5, 2, 1, 0), (5, 3, 2, 2), (5, 2, 2, 1),
    (6, 3, 3, 2), (6, 4, 4, 4), (6, 2, 2, 0),
    (7, 3, 2, 1), (7, 5, 5, 5), (7, 3, 3, 1),
    (8, 4, 3, 2), (8, 4, 4, 2), (8, 5, 1, 0),
    (9, 4, 3, 0), (9, 5, 3, 1), (9, 6, 6, 5),
    (10, 5, 5, 2), (11, 6, 4, 1),
]


def TestNumero2BosquesBiTerminalVsFuerzaBruta():
    """
        Test que compara, en 20 casos (n, a, b, c), la fórmula cerrada
        Numero2BosquesBiTerminal contra el cálculo por fuerza bruta
        NumeroDe2Bosques (menor del Laplaciano) sobre el grafo B(n; a, b, c).
    """
    if len(CASOS_BITERMINAL) != 20:
        return False
    for n, a, b, c in CASOS_BITERMINAL:
        if NumeroDe2Bosques(B(n, a, b, c)) != Numero2BosquesBiTerminal(n, a, b, c):
            return False
    return True


def TestEstructuraBiTerminal():
    """
        Test que verifica en los mismos casos que B(n; a, b, c) tiene
        n vértices, grados a y b en s y t, c vecinos comunes y el
        clique completo en el resto de los vértices.
    """
    for n, a, b, c in CASOS_BITERMINAL:
        TC = B(n, a, b, c)
        G = TC.grafo
        q = n - 2
        if G.order() != n or G.has_edge(TC.s, TC.t):
            return False
        if G.degree(TC.s) != a or G.degree(TC.t) != b:
            return False
        if len(set(G.neighbors(TC.s)) & set(G.neighbors(TC.t))) != c:
            return False
        if G.size() != q * (q - 1) // 2 + a + b:
            return False
    return True


def TestBiTerminalRechazaParametrosInvalidos():
    """
        Test que verifica que B lanza ValueError si no se cumple
        c <= b <= a o si a + b - c > n - 2.
    """
    for n, a, b, c in [(8, 2, 3, 1), (8, 3, 2, 3), (8, 3, 2, -1), (6, 4, 4, 1), (6, 5, 1, 0)]:
        try:
            B(n, a, b, c)
            return False
        except ValueError:
            pass
    return True


if __name__ == "__main__":
    tests = [
        ("TestNumero2BosquesBiTerminalVsFuerzaBruta", TestNumero2BosquesBiTerminalVsFuerzaBruta),
        ("TestEstructuraBiTerminal", TestEstructuraBiTerminal),
        ("TestBiTerminalRechazaParametrosInvalidos", TestBiTerminalRechazaParametrosInvalidos),
    ]

    fallidos = [nombre for nombre, test in tests if test() != True]
    if not fallidos:
        print("Todos los tests exitosos")
    else:
        for nombre in fallidos:
            print(f"Test {nombre} fallido")

import os
while not os.path.isdir("DataTypes"):
    os.chdir("..")

load("DataTypes/Grafo2Terminales.sage")
load("Scripts/GrafoBalloon.sage")
load("Scripts/VariantesGrafoBallon.sage")
load("Scripts/NumeroDe2Bosques.sage")
load("Scripts/Esqueleto.sage")
load("Scripts/DiferenciasT2.sage")

def TestDeltaT2HkVsFormula(k):
    """
        Test exhaustivo (n entre 4 y 20) que, para el índice k dado,
        compara Delta_k = t2(H_{n,m}^{(k)}) - t2(G_{n,m}) calculado
        numéricamente -- ambos números de 2-bosques vía el menor del
        Laplaciano, NumeroDe2Bosques -- contra evaluar directamente
        la fórmula cerrada de Delta_k (DeltaT2). Los (n, m) donde
        H_k no es construible, o donde (n', lambda') cae fuera del
        dominio probado de la fórmula, se saltean.
    """
    comparaciones = 0
    for n in range(4, 21):
        for m in range(n, binomial(n, 2) + 1):
            TC_original = ConstruirGrafo2TerminalesBalloon(n, m)
            b_valor, n_prima, lambda_prima = ParametrosEsqueletoBalloon(TC_original)

            try:
                formula = DeltaT2(k, b_valor, n_prima, lambda_prima)
            except ValueError:
                continue

            try:
                TC_Hk = VariantesBalloonH(k, n, m)
            except ValueError:
                continue

            numerico = NumeroDe2Bosques(TC_Hk) - NumeroDe2Bosques(TC_original)
            if numerico != formula:
                return False
            comparaciones += 1
    return comparaciones > 270


def TestDeltaT2H0VsFormula():
    """Delega en TestDeltaT2HkVsFormula(0), dedicado a k=0."""
    return TestDeltaT2HkVsFormula(0)


def TestDeltaT2H1VsFormula():
    """Delega en TestDeltaT2HkVsFormula(1), dedicado a k=1."""
    return TestDeltaT2HkVsFormula(1)


def TestDeltaT2H2VsFormula():
    """Delega en TestDeltaT2HkVsFormula(2), dedicado a k=2."""
    return TestDeltaT2HkVsFormula(2)


def TestDeltaT2RechazaFueraDeRango():
    """
        Test que verifica, con varios (n', lambda') fuera del
        dominio probado de cada fórmula, que DeltaT2H0/H1/H2 lanzan
        ValueError, y que DeltaT2 rechaza valores de k fuera de
        {0, 1, 2}.
    """
    casos = [
        (0, 4, 4), (0, 5, 5), (0, 5, 7), (0, 6, 6), (0, 6, 8),
        (1, 8, 8), (1, 4, 4), (1, 5, 5), (1, 6, 6), (1, 7, 7),
        (2, 4, 4), (2, 4, 6), (2, 5, 5), (2, 5, 7),
    ]
    b_valor = 5
    for k, n_prima, lambda_prima in casos:
        try:
            DeltaT2(k, b_valor, n_prima, lambda_prima)
            return False
        except ValueError:
            pass

    try:
        DeltaT2(3, b_valor, 10, 5)
        return False
    except ValueError:
        pass
    return True


if __name__ == "__main__":
    tests = [
        ("TestDeltaT2H0VsFormula", TestDeltaT2H0VsFormula),
        ("TestDeltaT2H1VsFormula", TestDeltaT2H1VsFormula),
        ("TestDeltaT2H2VsFormula", TestDeltaT2H2VsFormula),
        ("TestDeltaT2RechazaFueraDeRango", TestDeltaT2RechazaFueraDeRango),
    ]

    fallidos = [nombre for nombre, test in tests if test() != True]
    if not fallidos:
        print("Todos los tests exitosos")
    else:
        for nombre in fallidos:
            print(f"Test {nombre} fallido")

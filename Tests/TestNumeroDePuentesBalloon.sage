import os
while not os.path.isdir("DataTypes"):
    os.chdir("..")

load("DataTypes/Grafo2Terminales.sage")
load("Scripts/GrafoBalloon.sage")
load("Scripts/CantidadesEstructurales.sage")

def TestNumeroDePuentesCoincideConFormula():
    """
        Test exhaustivo (n entre 4 y 25) que compara, para todo
        (n, m) válido de un grafo Balloon G_{n,m}, la cantidad de
        puentes computada por fuerza bruta (b(TC), sobre el grafo
        con 2 terminales construido) contra la fórmula cerrada
        b(n, m).
    """
    comparaciones = 0
    for n in range(4, 26):
        for m in range(n, binomial(n, 2) + 1):
            TC = ConstruirGrafo2TerminalesBalloon(n, m)
            if b(TC) != b(n, m):
                return False
            comparaciones += 1
    return comparaciones > 400

def TestExtremos2PuentesImplicaDosPuentes():
    """
        Test exhaustivo (n entre 4 y 25) que verifica que, para
        todo m tal que primero <= m <= segundo (los extremos
        retornados por Extremos2Puentes(n)), el grafo Balloon
        G_{n,m} tiene exactamente 2 puentes.
    """
    probados = 0
    for n in range(4, 26):
        primero, segundo = Extremos2Puentes(n)
        for m in range(n, binomial(n, 2) + 1):
            if primero <= m <= segundo:
                probados += 1
                TC = ConstruirGrafo2TerminalesBalloon(n, m)
                if b(TC) != 2:
                    return False
    return probados > 100


if __name__ == "__main__":
    tests = [
        ("TestNumeroDePuentesCoincideConFormula", TestNumeroDePuentesCoincideConFormula),
        ("TestExtremos2PuentesImplicaDosPuentes", TestExtremos2PuentesImplicaDosPuentes),
    ]

    fallidos = [nombre for nombre, test in tests if test() != True]
    if not fallidos:
        print("Todos los tests exitosos")
    else:
        for nombre in fallidos:
            print(f"Test {nombre} fallido")

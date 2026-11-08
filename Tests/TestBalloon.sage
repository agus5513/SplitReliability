import os
while not os.path.isdir("DataTypes"):
    os.chdir("..")

load("DataTypes/Grafo2Terminales.sage")
load("Scripts/GrafoBalloon.sage")
load("Scripts/NumeroDe2Bosques.sage")

def TestNumeroDe2BosquesVsBallon():
    """
        Test exhaustivo (n entre 4 y 19), del mismo orden de
        magnitud que los tests de NumeroDe2BosquesH0/H1/H2 contra
        fuerza bruta, que compara, para todo grafo Balloon válido,
        el número de 2-bosques computado por fuerza bruta (vía el
        menor del Laplaciano) contra el obtenido por la fórmula
        exacta para grafos Balloon (NumeroDe2BosquesBallon).
    """

    comparaciones = 0
    for n in range(4, 20):
        for m in range(n, binomial(n, 2) + 1):
            TC = ConstruirGrafo2TerminalesBalloon(n, m)
            if NumeroDe2Bosques(TC) != NumeroDe2BosquesBallon(TC):
                return False
            comparaciones += 1
    return comparaciones > 400


if __name__ == "__main__":
    tests = [
        ("TestNumeroDe2BosquesVsBallon", TestNumeroDe2BosquesVsBallon),
    ]

    fallidos = [nombre for nombre, test in tests if test() != True]
    if not fallidos:
        print("Todos los tests exitosos")
    else:
        for nombre in fallidos:
            print(f"Test {nombre} fallido")

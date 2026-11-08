import os
while not os.path.isdir("DataTypes"):
    os.chdir("..")

load("DataTypes/Grafo2Terminales.sage")
load("Scripts/TerminalesSeparados.sage")
load("Scripts/GrafosRandom.sage")
load("Scripts/SplitReliability.sage")

import random

def TestSplitReliabilityVsArbol():
    """
        Test que compara el valor de la SplitReliability
        entre el implementación por fuerza bruta y la 
        fórmula de Brown. 
    """

    tol = 10**(-6)

    a = random.randint(10, 100)
    b = random.randint(10, 100)
    c = random.randint(10, 100)
    
    Arboles = [
        ConstruirArbol2TerminalesRandom(a),
        ConstruirArbol2TerminalesRandom(b),
        ConstruirArbol2TerminalesRandom(c)
    ]

    Probabilidades = [
        random.uniform(0, 1),
        random.uniform(0, 1),
        random.uniform(0, 1)
    ]
    
    for tree in Arboles: 
        for p in Probabilidades:
            if (SplitReliability(tree, p) - SplitReliabilityArbol(tree, p) > tol):
                return False
    return True
    

if __name__ == "__main__":
    if (TestSplitReliabilityVsArbol() == True):
        print("Test de SplitReliability contra árboles exitoso")
    else:
        print("Test de SplitReliability contra árboles fallido")




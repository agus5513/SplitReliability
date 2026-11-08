load("DataTypes/Grafo2Terminales.sage")
load("Scripts/TerminalesSeparados.sage")

from sage.all import Graph
from itertools import combinations

def N(TC : Grafo2Terminales, i):
    """
        Función que computa N_i(G): la cantidad de subgrafos
        generadores de G con i aristas que separan a los
        terminales s y t en dos componentes distintas.
    """
    G = TC.grafo
    m = G.size()

    if i < 0 or i > m:
        return 0

    acumulador = 0
    todas = set(G.edges())
    aristas = G.edges()
    for combo in combinations(aristas, i):
        a_mantener = set(combo)
        a_borrar = todas - a_mantener
        H = G.copy()
        H.delete_edges(a_borrar)
        AUX = Grafo2Terminales(H, TC.s, TC.t)
        if TerminalesSeparados(AUX):
            acumulador += 1
    return acumulador

def SplitReliability(TC : Grafo2Terminales, p):
    """
        Funcion que computa para un grafo de dos
        terminales, su SplitReliability para un
        parámetro p dado. Retorna un valor númerico,
        no el polinomio.
    """
    G = TC.grafo
    n = G.order()
    m = G.size()

    INICIO = n - 2
    FIN    = m - 1

    suma = 0
    for i in range(INICIO, FIN + 1):
        suma = suma + N(TC, i) * p**i * (1 - p)**(m-i)

    return suma

def SplitReliabilityArbol(TC : Grafo2Terminales, p):
    """
        Función que computa la SplitReliability
        numérica de un árbol.
    """
    k = TC.grafo.distance(TC.s, TC.t)
    n = TC.grafo.order()
    return k * p**(n - 2) * (1 - p)

def SplitReliabilityCiclo(TC : Grafo2Terminales, p):
    """
        Función que computa la SplitReliability
        numérica de un grafo cíclico.
    """
    k = TC.grafo.distance(TC.s, TC.t)
    n = TC.grafo.order()
    return k * (n - k) * p**(n - 2) * (1 - p)**2

def F(TC : Grafo2Terminales, i):
    """
        Función que computa F_i(G) = N_{m-i}(G), según
        la definición introducida en RomeroTCS-2025.

        F_i(G) cuenta la cantidad de subgrafos generadores
        de G con m - i aristas que separan a los terminales
        s y t en dos componentes distintas (equivalentemente,
        subgrafos "split" obtenidos al eliminar i aristas).
    """
    m = TC.grafo.size()
    j = m - i
    return N(TC, j)

def VectorF(TC : Grafo2Terminales):
    """
        Función que computa el vector F(G) = (F_0(G), F_1(G), ..., F_m(G)),
        según RomeroTCS-2025.
    """
    G = TC.grafo
    m = G.size()
    return [F(TC, i) for i in range(0, m + 1)]









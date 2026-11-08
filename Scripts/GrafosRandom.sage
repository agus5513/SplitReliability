load("DataTypes/Grafo2Terminales.sage")

from sage.all import Graph
import random

def SimularCaida(G : Grafo2Terminales, p):
    """
        Función que recibe el estructurado de un grafo de
        2 terminales y simula que cada arista se mantiene
        operativa con una probabilidad de p (o equivalentemete
        que cada arista se elimina con probabilidad 1-p)
        Retorna un nuevo estructurado sin afectar el parámetro
    """
    H = G.grafo.copy()
    aristas = H.edges(labels=False)
    for e in aristas:
        if (random.uniform(0, 1) > p):
            H.delete_edge(e)
    return Grafo2Terminales(H, G.s, G.t)

def ConstruirArbol2TerminalesRandom(n):
    T = graphs.RandomTree(n)
    s = 0
    t = random.randint(1, n - 1)
    return Grafo2Terminales(T, s, t)

def ConstruirCiclo2TerminalesRandom(n):
    C = graphs.CycleGraph(n)
    s = 0
    t = random.randint(1, n - 1)
    return Grafo2Terminales(C, s, t)

def ConstruirGrafo2TerminalesRandomSinAristaST(n, m, intentos=1000):
    """
        Función que genera un grafo aleatorio (Erdos-Renyi
        G(n,m)) conexo, con terminales s=0, t=n-1, tal que
        (s,t) no es una arista del grafo. Reintenta hasta
        `intentos` veces antes de fallar.
    """
    s, t = 0, n - 1
    for _ in range(intentos):
        G = graphs.RandomGNM(n, m)
        if G.is_connected() and not G.has_edge(s, t):
            return Grafo2Terminales(G, s, t)
    raise RuntimeError("No se pudo generar un grafo aleatorio conexo sin la arista (s,t)")
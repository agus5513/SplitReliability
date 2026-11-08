load("DataTypes/Grafo2Terminales.sage")

from sage.all import Graph

def ParticionMultipartito(p, e):
    """
        Función para uso interno. Dada la cantidad de vértices p
        y el exceso total e (en {-2, -1, 0}), retorna una lista
        decreciente de partes en {1,2,3,4} que suman p y tiene
        exceso total e; el resto se rellena con partes de tamaño 3.
        Retorna None si no existe tal partición (solo ocurre
        para e = -1, p = 6, es decir el caso n = 7).
    """
    ajustes = {-2: {0: [2, 1], 1: [2, 2], 2: [1, 1]},
               -1: {0: [4, 2, 2, 1], 1: [1], 2: [2]},
                0: {0: [], 1: [4, 2, 1], 2: [4, 2, 2]}}
    a = ajustes[e][p % 3]
    resto = p - sum(a)
    if resto < 0:
        return None
    return sorted(a + [3] * (resto // 3), reverse=True)

def ConstruirGrafoMejorEn2Puentes(n, m):
    """
        Función que construye el grafo con 2 terminales de n vértices
        y m aristas que se propone como mejor candidato en el caso de
        2 puentes. Consiste en un grafo multipartito completo con
        p = n - 1 vértices al que se le agrega un vértice t unido por
        una única arista a un vértice x de la parte más grande.
        El terminal s se ubica en la misma parte que x (x != s).
    """
    n1, m1 = n - 2, m - 2
    lam = m1 - binomial(n1 - 1, 2)
    kp = n1 - 1 - lam
    p = n - 1

    if kp - 2 not in (-2, -1, 0):
        raise ValueError(f"La construcción no está definida para n={n}, m={m}")

    partes = ParticionMultipartito(p, kp - 2)
    if partes is None:
        raise ValueError(f"No existe la construcción para n={n}, m={m}")

    G = graphs.CompleteMultipartiteGraph(partes)  # vértices 0..p-1, parte 0 = la más grande
    x, s, t = 0, 1, p
    G.add_vertex(t)
    G.add_edge(x, t)

    if G.size() != m:
        raise ValueError(f"La construcción tiene {G.size()} aristas y se esperaban {m}")
    return Grafo2Terminales(G, s, t)

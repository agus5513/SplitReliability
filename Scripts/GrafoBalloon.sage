load("DataTypes/Grafo2Terminales.sage")

from sage.all import Graph

def Indice(n, m):
    return n >= 4 and n <= m <= binomial(n, 2)

def Indice0(n, m):
    return n >= 4 and binomial(n - 1, 2) + 2 <= m <= binomial(n, 2)

def Indice1(n, m):
    return Indice(n, m) and not Indice0(n, m)

def EstructuraBalloon(n, m):
    """
        Función recursiva que construye la estructura
        del grafo Balloon. Función para uso interno en
        la librería ya que no coloca terminales s,t.
    """
    if Indice0(n, m):
        G = graphs.CompleteGraph(n - 1)
        extremos = m - binomial(n - 1, 2)
        for v in range(extremos):
            G.add_edge(n - 1, v)
    elif n == 4 and m == 4:
        G = graphs.CompleteGraph(3)
        G.add_edge(n - 1, 0)
    else:
        G = EstructuraBalloon(n - 1, m - 1)
        v_min = min(G.vertices(), key=G.degree)
        G.add_edge(n - 1, v_min)
    return G

def ConstruirGrafo2TerminalesBalloon(n, m):
    """
        Función que computa el grafo Balloon con 2
        terminales para un valor de n,m.
    """
    G = EstructuraBalloon(n, m)
    s, t = 0, 0
    distancia_maxima = 0
    for u in G.vertices():
        for v in G.vertices():
            d = G.distance(u, v)
            if d > distancia_maxima:
                s, t = u, v
                distancia_maxima = d
    return Grafo2Terminales(G, s, t)

load("DataTypes/Grafo2Terminales.sage")
load("Scripts/CantidadesEstructurales.sage")

from sage.all import Graph

def ObtenerEsqueleto(TC):
    """
    Recibe como parámetro un grafo con 2 terminales
    y retorna un nuevo grafo con 2 terminales que se obtiene
    tras colapsar todos los puentes del grafo original. Proyecta
    los vértices terminales.
    """

    G = TC.grafo
    puentes = G.bridges(labels=False)

    padre = {v: v for v in G.vertices()}

    def encontrar(v):
        while padre[v] != v:
            padre[v] = padre[padre[v]]
            v = padre[v]
        return v

    def unir(u, v):
        ru, rv = encontrar(u), encontrar(v)
        if ru != rv:
            padre[ru] = rv

    for u, v in puentes:
        unir(u, v)

    H = Graph()
    H.add_vertices({encontrar(v) for v in G.vertices()})
    for u, v in G.edges(labels=False):
        ru, rv = encontrar(u), encontrar(v)
        if ru != rv:
            H.add_edge(ru, rv)

    s_proyectado = encontrar(TC.s)
    t_proyectado = encontrar(TC.t)
    return Grafo2Terminales(H, s_proyectado, t_proyectado)

def ParametrosEsqueletoBalloon(TC):
    """
        Función auxiliar que retorna (b, n', lambda') para el grafo
        Balloon TC: la cantidad de puentes de TC y el orden y grado
        mínimo del esqueleto de TC (mismos parámetros que usa
        NumeroDe2BosquesVarianteBallon).
    """
    b_valor = b(TC)
    esqueleto = ObtenerEsqueleto(TC)
    G_prima = esqueleto.grafo
    n_prima = G_prima.order()
    m_prima = G_prima.size()
    lambda_prima = m_prima - binomial(n_prima - 1, 2)
    return b_valor, n_prima, lambda_prima

load("DataTypes/Grafo2Terminales.sage")

from sage.all import Graph

def TerminalesOptimos(G: Graph, p):
    """
        Función que recibe como parámetros una
        probabilidad p y un grafo G, que computa
        las ubicaciones de s, t óptimas para
        maximizar su SplitReliability.
    """

    n = G.order()
    MaxSplitReliability = -1
    UbicacionOptima = None

    for s_opt in range(0, n):
        for t_opt in range(0, n):
            Candidato = Grafo2Terminales(G, s_opt, t_opt)
            SR = SplitReliability(Candidato, p)
            if (SR > MaxSplitReliability):
                UbicacionOptima = Grafo2Terminales(G, s_opt, t_opt)
                MaxSplitReliability = SR
    return UbicacionOptima
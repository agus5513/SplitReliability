load("DataTypes/Grafo2Terminales.sage")

from sage.all import Graph

def TerminalesSeparados(G : Grafo2Terminales):
    """
        Función booleana que determina si el grafo cuenta
        con exactamente 2 componenetes conexas con un único
        nodo terminal en cada una de ellas
    """
    # Computamos las componenes conexas del grafo 
    comp = G.grafo.connected_components()
    if len(comp) != 2:
        return False
    if (G.s in comp[0] and G.t in comp[1] or G.s in comp[1] and G.t in comp[0]):
        return True
    return False
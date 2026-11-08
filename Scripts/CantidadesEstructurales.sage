from sage.all import Graph

def b(*args):
    """
        b(G): computa la cantidad de puentes (bridges) de G,
        es decir, aristas cuya eliminación aumenta la cantidad
        de componentes conexas de G. G es un Grafo2Terminales.

        b(n, m): computa la cantidad de puentes del grafo
        Balloon G_{n,m} (n vértices, m aristas) vía la fórmula
        cerrada n - 3 - floor((sqrt(8*(m - n) + 1) - 1) / 2),
        sin necesidad de construir el grafo. Coincide con
        b(G) para el grafo Balloon correspondiente, y con el
        máximo de puentes entre todos los grafos conexos de
        n vértices y m aristas (verificado por fuerza bruta
        contra GrafosQueMaximizanB).
    """
    if len(args) == 1:
        G, = args
        return len(list(G.grafo.bridges()))
    if len(args) == 2:
        n, m = args
        t = floor((sqrt(8*(m - n) + 1) - 1) / 2)
        return n - 3 - t
    raise TypeError("b espera 1 argumento (G) o 2 argumentos (n, m)")

def Extremos2PuentesRestriccionLambda(n):
    """
        Función que retorna, para un grafo Balloon de n
        vértices, los extremos (primero, segundo) del intervalo
        de valores de m tales que primero <= m <= segundo implica
        que G_{n,m} tiene exactamente 2 puentes.
    """
    primero = n + (n - 5)*(n - 4)/2
    segundo = n + (n - 4)*(n - 3)/2 - 4

    return primero, segundo

def Lambda(G):
    """
        Función que computa lambda(G): la mínima cantidad
        de aristas que deben retirarse de G para aumentar
        su cantidad de componentes conexas (conectividad
        de aristas). Wrapper directo de edge_connectivity.
    """
    return G.grafo.edge_connectivity()

def t(G):
    """
        Función que computa t(G): el número de árboles
        recubridores de G. Wrapper directo de
        number_of_spanning_trees.
    """
    return G.grafo.spanning_trees_count()

def GenerarGrafosConexos(n, m):
    """
        Función auxiliar que genera, vía Nauty, todos los
        grafos conexos no isomorfos de n vértices y m aristas.
    """
    return graphs.nauty_geng(f"{n} -c {m}:{m}")

def GrafosQueMaximizanB(n, m):
    """
        Función que busca, entre todos los grafos conexos
        de n vértices y m aristas, aquellos que maximizan
        b(G). Retorna una tupla (valor_maximo, lista_de_grafos).
    """
    valor_maximo = -1
    optimos = []
    for G in GenerarGrafosConexos(n, m):
        valor = b(G)
        if valor > valor_maximo:
            valor_maximo = valor
            optimos = [G]
        elif valor == valor_maximo:
            optimos.append(G)
    return valor_maximo, optimos

def GrafosQueMinimizanLambda(n, m):
    """
        Función que busca, entre todos los grafos conexos
        de n vértices y m aristas, aquellos que minimizan
        lambda(G). Retorna una tupla (valor_minimo, lista_de_grafos).
    """
    valor_minimo = None
    optimos = []
    for G in GenerarGrafosConexos(n, m):
        valor = Lambda(G)
        if valor_minimo is None or valor < valor_minimo:
            valor_minimo = valor
            optimos = [G]
        elif valor == valor_minimo:
            optimos.append(G)
    return valor_minimo, optimos

def GrafosQueMinimizanT(n, m):
    """
        Función que busca, entre todos los grafos conexos
        de n vértices y m aristas, aquellos que minimizan
        t(G). Retorna una tupla (valor_minimo, lista_de_grafos).
    """
    valor_minimo = None
    optimos = []
    for G in GenerarGrafosConexos(n, m):
        valor = t(G)
        if valor_minimo is None or valor < valor_minimo:
            valor_minimo = valor
            optimos = [G]
        elif valor == valor_minimo:
            optimos.append(G)
    return valor_minimo, optimos

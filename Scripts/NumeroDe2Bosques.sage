load("DataTypes/Grafo2Terminales.sage")
load("Scripts/Esqueleto.sage")
load("Scripts/CantidadesEstructurales.sage")

from sage.all import Graph

def NumeroDe2Bosques(TC):
    """
    Cuenta la cantidad de bosques de tamaño 2 (arboles disjuntos)
    con s y t en árboles distintos, vía el menor del Laplaciano.
    """

    G = TC.grafo
    L = G.laplacian_matrix()
    n = G.vertices()
    idx = [i for i, v in enumerate(n) if v != TC.s and v != TC.t]
    M = L.matrix_from_rows_and_columns(idx, idx)
    return M.determinant()

def NumeroDe2BosquesBallon(TC):
    """
    Utiliza la fórmula exacta para el número de 2-bosques
    pero es solo para grafos con estructura balloon...
    """

    esqueleto = ObtenerEsqueleto(TC)
    G_prima = esqueleto.grafo
    n_prima = G_prima.order()
    m_prima = G_prima.size()

    lambda_prima = m_prima - binomial(n_prima - 1, 2)
    r_prima = n_prima - 1 - lambda_prima

    arboles_esqueleto = lambda_prima * (n_prima)**(lambda_prima - 1) * (n_prima - 1)**(n_prima - lambda_prima - 2)
    
    if r_prima > 0:
        formula_esqueleto = (2*n_prima - r_prima - 1) * (n_prima)**(n_prima - r_prima - 2) * (n_prima - 1)**(r_prima - 2)
    else:
        formula_esqueleto = 2 * (n_prima)**(n_prima - 3)

    return b(TC) * arboles_esqueleto + formula_esqueleto

def GrafosQueMaximizanT2(n, m):
    """
        Función que busca, por fuerza bruta entre todos los grafos
        conexos de n vértices y m aristas (vía GenerarGrafosConexos)
        y todas las ubicaciones de terminales s,t, aquellos
        Grafo2Terminales que maximizan T2 (NumeroDe2Bosques).
        Retorna una tupla (valor_maximo, lista_de_grafos).
    """
    valor_maximo = -1
    optimos = []
    for G in GenerarGrafosConexos(n, m):
        vertices = G.vertices()
        for i in range(len(vertices)):
            for j in range(i + 1, len(vertices)):
                TC = Grafo2Terminales(G, vertices[i], vertices[j])
                valor = NumeroDe2Bosques(TC)
                if valor > valor_maximo:
                    valor_maximo = valor
                    optimos = [TC]
                elif valor == valor_maximo:
                    optimos.append(TC)
    return valor_maximo, optimos
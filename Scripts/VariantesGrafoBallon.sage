load("DataTypes/Grafo2Terminales.sage")
load("Scripts/GrafoBalloon.sage")
load("Scripts/Esqueleto.sage")
load("Scripts/CantidadesEstructurales.sage")

from sage.all import Graph

def RepresentantesEsqueleto(TC : Grafo2Terminales):
    """
        Función auxiliar que computa, para cada vértice de TC.grafo,
        su representante en el esqueleto G' tras colapsar los
        puentes (mismo criterio de unión utilizado en
        ObtenerEsqueleto). Permite recuperar la arista real de G
        asociada a una arista de G'.
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

    return {v: encontrar(v) for v in G.vertices()}

def NodoSPrima(TC : Grafo2Terminales):
    """
        Función auxiliar que retorna una tupla (G', s', lambda')
        donde G' es el esqueleto de TC, s' es el vértice de grado
        mínimo de G' y lambda' es el grado de dicho vértice.
    """
    esqueleto = ObtenerEsqueleto(TC)
    G_prima = esqueleto.grafo
    s_prima = min(G_prima.vertices(), key=G_prima.degree)
    lambda_prima = G_prima.degree(s_prima)
    return G_prima, s_prima, lambda_prima

def AristaRealDeEsqueleto(TC : Grafo2Terminales, representantes, a, b):
    """
        Función auxiliar que retorna una arista real de TC.grafo
        cuya proyección en el esqueleto es la arista (a, b) de G'.
    """
    for u, v in TC.grafo.edges(labels=False):
        if {representantes[u], representantes[v]} == {a, b}:
            return (u, v)
    return None

def ContraerPuente(TC : Grafo2Terminales):
    """
        Función auxiliar que contrae un puente de TC.grafo (el
        primero reportado por Sage) y retorna un nuevo
        Grafo2Terminales, actualizando los terminales en caso de
        que alguno de ellos haya sido fusionado por la contracción.
    """
    G = TC.grafo.copy()
    puentes = list(G.bridges(labels=False))
    if not puentes:
        raise ValueError("G no tiene puentes para contraer")

    u, v = puentes[0]
    G.contract_edge((u, v))
    s = u if TC.s == v else TC.s
    t = u if TC.t == v else TC.t
    return Grafo2Terminales(G, s, t)

def H0(n, m):
    TC = ContraerPuente(ConstruirGrafo2TerminalesBalloon(n, m))
    G_prima, s_prima, lambda_prima = NodoSPrima(TC)

    representantes = RepresentantesEsqueleto(TC)
    vecino = next(iter(G_prima.neighbors(s_prima)))
    u, v = AristaRealDeEsqueleto(TC, representantes, s_prima, vecino)

    G = TC.grafo.copy()
    G.subdivide_edge(u, v, 1)
    return Grafo2Terminales(G, TC.s, TC.t)

def H1(n, m):
    TC = ContraerPuente(ConstruirGrafo2TerminalesBalloon(n, m))
    G_prima, s_prima, lambda_prima = NodoSPrima(TC)

    representantes = RepresentantesEsqueleto(TC)
    s_real_prima = representantes[TC.s]

    no_adyacentes = [
        v for v in G_prima.vertices()
        if v != s_prima and v != s_real_prima and not G_prima.has_edge(s_prima, v)
    ]
    if len(no_adyacentes) < 2:
        raise ValueError("G' no tiene 2 nodos no adyacentes a s' (excluyendo la imagen de s)")

    a, b = no_adyacentes[0], no_adyacentes[1]
    if not G_prima.has_edge(a, b):
        raise ValueError("los nodos elegidos no están unidos por una arista en G'")

    representantes_edge = AristaRealDeEsqueleto(TC, representantes, a, b)
    u, v = representantes_edge

    G = TC.grafo.copy()
    G.subdivide_edge(u, v, 1)
    return Grafo2Terminales(G, TC.s, TC.t)

def H2(n, m):
    TC = ContraerPuente(ConstruirGrafo2TerminalesBalloon(n, m))
    G_prima, s_prima, lambda_prima = NodoSPrima(TC)

    adyacentes = [v for v in G_prima.neighbors(s_prima) if v != s_prima]
    if len(adyacentes) < 2:
        raise ValueError("s' no tiene 2 vecinos en G'")

    a, b = adyacentes[0], adyacentes[1]
    if not G_prima.has_edge(a, b):
        raise ValueError("los vecinos de s' elegidos no están unidos por una arista en G'")

    representantes = RepresentantesEsqueleto(TC)
    u, v = AristaRealDeEsqueleto(TC, representantes, a, b)

    G = TC.grafo.copy()
    G.subdivide_edge(u, v, 1)
    return Grafo2Terminales(G, TC.s, TC.t)

def VariantesBalloonH(k, n, m):
    """
        Función wrapper que retorna el grafo H_k(n, m) construido
        según H0, H1 o H2, de acuerdo al valor de k en {0, 1, 2}.
    """
    if k == 0:
        return H0(n, m)
    elif k == 1:
        return H1(n, m)
    elif k == 2:
        return H2(n, m)
    raise ValueError(f"k={k} no es un índice válido (debe ser 0, 1 o 2)")

def NumeroDe2BosquesH0(b_valor, n_prima, lambda_prima):
    """
        Fórmula exacta para N_{n-2}(H_{n,m}^{(0)}), válida en el
        dominio 2 < lambda' <= n' - 2.
    """
    if not (2 < lambda_prima <= n_prima - 2):
        raise ValueError(
            f"NumeroDe2BosquesH0: fórmula no definida fuera de "
            f"2 < lambda' <= n' - 2 (se obtuvo n'={n_prima}, lambda'={lambda_prima})"
        )
    return (
        (b_valor - 1) * (
            lambda_prima * n_prima**(lambda_prima - 1) * (n_prima - 1)**(n_prima - lambda_prima - 2)
            + (lambda_prima - 1) * n_prima**(lambda_prima - 2) * (n_prima - 1)**(n_prima - lambda_prima - 1)
        )
        + (
            (n_prima + lambda_prima) * n_prima**(lambda_prima - 1) * (n_prima - 1)**(n_prima - lambda_prima - 3)
            + (n_prima + lambda_prima - 1) * n_prima**(lambda_prima - 2) * (n_prima - 1)**(n_prima - lambda_prima - 2)
        )
    )

def NumeroDe2BosquesH1(b_valor, n_prima, lambda_prima):
    """
        Fórmula exacta para N_{n-2}(H_{n,m}^{(1)}), válida en el
        dominio lambda' <= n' - 4.
    """
    if lambda_prima > n_prima - 4:
        raise ValueError(
            f"NumeroDe2BosquesH1: fórmula no definida "
            f"para lambda'={lambda_prima} > n' - 4={n_prima - 4}"
        )
    return (
        (b_valor - 1) * (
            lambda_prima * n_prima**(lambda_prima - 1) * (n_prima - 1)**(n_prima - lambda_prima - 2)
            + lambda_prima * n_prima**(lambda_prima - 1) * (n_prima - 3) * (n_prima - 1)**(n_prima - lambda_prima - 3)
        )
        + (
            (n_prima + lambda_prima) * n_prima**(lambda_prima - 1) * (n_prima - 1)**(n_prima - lambda_prima - 3)
            + (n_prima - 3) * (n_prima + lambda_prima) * n_prima**(lambda_prima - 1) * (n_prima - 1)**(n_prima - lambda_prima - 4)
        )
    )

def NumeroDe2BosquesH2(b_valor, n_prima, lambda_prima):
    """
        Fórmula exacta para N_{n-2}(H_{n,m}^{(2)}), válida en el
        dominio 2 <= lambda' <= n' - 2.
    """
    if not (2 <= lambda_prima <= n_prima - 2):
        raise ValueError(
            f"NumeroDe2BosquesH2: fórmula no definida fuera de "
            f"2 <= lambda' <= n' - 2 (se obtuvo n'={n_prima}, lambda'={lambda_prima})"
        )
    return (
        (b_valor - 1) * 2 * lambda_prima * n_prima**(lambda_prima - 2) * (n_prima - 1)**(n_prima - lambda_prima - 1)
        + 2 * (n_prima + lambda_prima) * n_prima**(lambda_prima - 2) * (n_prima - 1)**(n_prima - lambda_prima - 2)
    )

def NumeroDe2BosquesVarianteBallon(k, n, m):
    """
        Función wrapper que computa, vía fórmula exacta, el número
        de 2-bosques de H_{n,m}^{(k)} (k = 0, 1, 2), sin construir
        el grafo. Delega en NumeroDe2BosquesH0, NumeroDe2BosquesH1
        o NumeroDe2BosquesH2 según el valor de k, pasándoles
        b = b(G_{n,m}) (cantidad de puentes del grafo Balloon
        G_{n,m}), n' (orden del esqueleto G' de G_{n,m}) y lambda'
        (grado mínimo de un nodo de G'). Cada una de esas funciones
        lanza ValueError si (n', lambda') queda fuera de su dominio
        de validez.
    """
    TC = ConstruirGrafo2TerminalesBalloon(n, m)
    b_valor = b(TC)

    esqueleto = ObtenerEsqueleto(TC)
    G_prima = esqueleto.grafo
    n_prima = G_prima.order()
    m_prima = G_prima.size()
    lambda_prima = m_prima - binomial(n_prima - 1, 2)

    if k == 0:
        return NumeroDe2BosquesH0(b_valor, n_prima, lambda_prima)
    elif k == 1:
        return NumeroDe2BosquesH1(b_valor, n_prima, lambda_prima)
    elif k == 2:
        return NumeroDe2BosquesH2(b_valor, n_prima, lambda_prima)
    raise ValueError(f"k={k} no es un índice válido (debe ser 0, 1 o 2)")

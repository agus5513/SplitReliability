r"""
TeoremaK3.sage

Modelos G*_{n,m} del Teorema A (k = 3, b(n,m) >= 2), formula de descomposicion
y verificacion contra la busqueda exhaustiva.

Requiere DataTypes/GrafoKTerminales.sage y Scripts/FiltroLexicografico.sage.
"""

from itertools import combinations


def b_nm(n, m):
    r"""b(n,m) para m >= n (formula corregida); para m = n - 1 vale n - 1."""
    return maximo_puentes_formula(n, m)


def terna_justa(b):
    r"""(l1, l2, l3) justa con l1 >= l2 >= l3 y suma b."""
    q, r = divmod(b, 3)
    return tuple(sorted([q + (1 if i < r else 0) for i in range(3)], reverse=True))


def nucleo_globo(np_, mp_):
    r"""
    Globo B_{n',m'} casi completo: K_{n'-1} en 0..n'-2 y z = n'-1 adyacente a
    0..lambda'-1. Devuelve (B, z, lambda').
    """
    lam = mp_ - binomial(np_ - 1, 2)
    B = graphs.CompleteGraph(np_ - 1)
    B.add_vertex(np_ - 1)
    for u in range(lam):
        B.add_edge(np_ - 1, u)
    return B, np_ - 1, lam


def theta(H, patas):
    r"""
    Theta(H; (u1,l1), (u2,l2), (u3,l3)): a cada u_r se le cuelga un camino de
    largo l_r terminado en un terminal nuevo (si l_r = 0 el terminal es u_r).
    """
    G = H.copy(immutable=False)
    sig = max(G.vertices(sort=True)) + 1
    terminales = []
    for u, l in patas:
        actual = u
        for _ in range(l):
            G.add_edge(actual, sig)
            actual = sig
            sig += 1
        terminales.append(actual)
    return GrafoKTerminales(G, terminales)


def modelo(n, m):
    r"""G*_{n,m} del Teorema A (requiere b(n,m) >= 1 para tener sentido como plantilla)."""
    b = b_nm(n, m)
    np_, mp_ = n - b, m - b
    B, z, lam = nucleo_globo(np_, mp_)
    l1, l2, l3 = terna_justa(b)
    if lam == np_ - 1:
        return theta(B, [(0, l1), (1, l2), (2, l3)])
    if lam == np_ - 2:
        p, q = z, np_ - 2          # arista faltante pq
        if l2 == l3:
            return theta(B, [(0, l1), (p, l2), (q, l3)])
        return theta(B, [(p, l1), (p, l2), (q, l3)])
    w1, w2 = lam, lam + 1          # dos no vecinos de z
    return theta(B, [(w1, l1), (w2, l2), (z, l3)])


def variantes_clase(n, m):
    r"""
    En el caso lambda' = n'-2 con l2 > l3 (b = 2 mod 3, b >= 5): miembros no
    isomorfos de la clase de G*: nucleo en la hoja de la pata corta y, si l3 >= 2,
    nucleo en el interior de la pata corta.
    """
    b = b_nm(n, m)
    np_, mp_ = n - b, m - b
    B, z, lam = nucleo_globo(np_, mp_)
    l1, l2, l3 = terna_justa(b)
    if not (lam == np_ - 2 and l2 > l3 and b >= 3):
        return []
    p, q = z, np_ - 2
    salida = []
    for d in range(1, l3 + 1):     # d = distancia del nucleo a la mediana sobre la pata corta
        G = B.copy(immutable=False)
        sig = max(G.vertices(sort=True)) + 1
        # camino desde q hasta la mediana v, de largo d
        actual = q
        for _ in range(d):
            G.add_edge(actual, sig)
            actual = sig
            sig += 1
        v = actual
        ts = []
        for l in (l1, l2):
            actual = v
            for _ in range(l):
                G.add_edge(actual, sig)
                actual = sig
                sig += 1
            ts.append(actual)
        # resto de la pata corta: desde p, largo l3 - d (si d = l3, el terminal es p)
        actual = p
        for _ in range(l3 - d):
            G.add_edge(actual, sig)
            actual = sig
            sig += 1
        ts.append(actual)
        salida.append(GrafoKTerminales(G, ts))
    return salida


def verificar(n, m, verbose=True):
    r = ganadores(n, m, 3)
    G = modelo(n, m)
    ok = G.F_tupla() == r.tupla_maxima
    clase_ok = all(H.F_tupla() == r.tupla_maxima for H in variantes_clase(n, m))
    iso = [H.es_isomorfo(G) for H in r.ganadores]
    if verbose:
        b = b_nm(n, m)
        np_ = n - b
        lam = m - b - binomial(np_ - 1, 2)
        print(f"({n},{m}) b={b} n'={np_} lam'={lam} legs={terna_justa(b)}: "
              f"modelo {'OK' if ok else 'FALLA'}, variantes {'OK' if clase_ok else 'FALLA'}, "
              f"clases iso ganadoras={len(r.ganadores)}, decide en F_{r.indice_discriminante}",
              flush=True)
    return ok and clase_ok, r


# ----------------------------------------------------------------------
# Formula de descomposicion (Proposicion de factorizacion)
# ----------------------------------------------------------------------

def coef_conexos(C, j):
    r"""c_j(C): # conjuntos de j aristas cuya remocion deja C conexo."""
    E = list(C.edges(labels=False))
    total = 0
    for X in combinations(E, j):
        H = C.copy(immutable=False)
        H.delete_edges(X)
        total += H.is_connected()
    return total


def formula_descomposicion(H):
    r"""
    Para H en T(2) (b(H) >= 2, un unico nucleo), calcula F_i(H) por la formula
    phi_b c_{i-2}(C) + sum_e s_{i-1}(C; u_e, v_e) + F_i(C; pi(t)).
    """
    G = H.grafo()
    puentes = H.puentes()
    b = len(puentes)
    nucleos = H.nucleos()
    assert len(nucleos) == 1
    C = nucleos[0]
    VC = set(C.vertices(sort=False))
    m = G.size()
    T = G.copy(immutable=False)
    T.delete_edges(C.edges(labels=False))

    def proy(v):
        comp = T.connected_component_containing_vertex(v)
        return next(x for x in comp if x in VC)

    pi = [proy(t) for t in H.terminales()]
    ell = sorted([0, 0, 0])
    l = terna_justa(b)
    phi = l[0] * l[1] + l[0] * l[2] + l[1] * l[2]
    F = [0] * (m + 1)
    Ctab = TablaConexos(C)
    # termino con dos puentes
    mC = C.size()
    for j in range(mC + 1):
        if j + 2 <= m:
            F[j + 2] += phi * coef_conexos(C, j)
    # termino con un puente
    for e in puentes:
        G2 = G.copy(immutable=False)
        G2.delete_edge(e)
        lado = [c for c in G2.connected_components() if not (set(c) & VC)][0]
        dentro = [t for t in H.terminales() if t in lado]
        if len(dentro) != 1:
            continue
        otros = [pi[i] for i, t in enumerate(H.terminales()) if t not in lado]
        if otros[0] == otros[1]:
            continue
        s = GrafoKTerminales(C, otros, tabla=Ctab).F_polinomio()
        for j, c in enumerate(s):
            if c:
                F[j + 1] += c
    # termino sin puentes
    if len(set(pi)) == 3:
        f3 = GrafoKTerminales(C, pi, tabla=Ctab).F_polinomio()
        for j, c in enumerate(f3):
            F[j] += c
    return F

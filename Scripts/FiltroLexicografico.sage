r"""
FiltroLexicografico.sage

Busqueda exhaustiva de los grafos k-terminales localmente mas separables
(F-tupla lexicograficamente maxima) y verificaciones de referencia:

  - verificar_teorema10(n_max): reproduce el Teorema 10 de Romero (2025),
    k = 2: el globo G_{n,m} tiene la F-tupla maxima en cada (n, m) de I.
  - tabla_k3(ns): ganadores para k = 3, comparados con los prefijos obtenidos
    en la exploracion inicial (n = 6, 7).
  - verificar_menor_arbol(n_max, ks): F_{k-1}(T) = det L(T)_[R] en arboles.

Requiere haber cargado antes DataTypes/GrafoKTerminales.sage.

Uso:
    load("DataTypes/GrafoKTerminales.sage")
    load("Scripts/FiltroLexicografico.sage")
    r = ganadores(7, 10, 3)
    r.resumen()
"""

from itertools import combinations

try:
    GrafoKTerminales
except NameError:
    raise ImportError("Carga primero DataTypes/GrafoKTerminales.sage")


# ----------------------------------------------------------------------
# Enumeracion
# ----------------------------------------------------------------------

def grafos_conexos(n, m):
    r"""Todos los grafos simples conexos con n vertices y m aristas, sin repetir isomorfos (nauty)."""
    if n == 1:
        if m == 0:
            yield Graph(1)
        return
    yield from graphs.nauty_geng(f"{n} {m}:{m} -c")


def candidatos(G, k):
    r"""
    Un GrafoKTerminales por cada orbita de k-subconjuntos de vertices de G
    bajo Aut(G). Todos comparten la misma TablaConexos.
    """
    tabla = TablaConexos(G)
    vistos = set()
    for terminales in combinations(G.vertices(sort=True), k):
        H = GrafoKTerminales(G, terminales, tabla=tabla)
        clave = H.forma_canonica()
        if clave not in vistos:
            vistos.add(clave)
            yield H


class ResultadoFiltro:
    r"""
    Resultado de ganadores(n, m, k).

    Atributos:
        - n, m, k
        - tupla_maxima: F-tupla (F_1, ..., F_m) lexicograficamente maxima.
        - ganadores: lista de GrafoKTerminales (uno por clase de isomorfismo)
          con esa F-tupla. Son todos split-equivalentes entre si.
        - segunda_tupla: la mejor F-tupla distinta de la maxima (o None).
        - indice_discriminante: menor i con F_i(ganador) != F_i(segundo);
          es el indice donde se decide la clase ganadora.
        - maximo_puentes: b(n, m) contado directamente sobre la clase.
        - candidatos: cuantos grafos k-terminales no isomorfos se compararon.
    """

    def __init__(self, n, m, k, tupla_maxima, ganadores, segunda_tupla,
                 maximo_puentes, candidatos):
        self.n, self.m, self.k = n, m, k
        self.tupla_maxima = tupla_maxima
        self.ganadores = ganadores
        self.segunda_tupla = segunda_tupla
        self.maximo_puentes = maximo_puentes
        self.candidatos = candidatos
        self.indice_discriminante = None
        if segunda_tupla is not None:
            self.indice_discriminante = next(
                i + 1 for i, (a, b) in enumerate(zip(tupla_maxima, segunda_tupla))
                if a != b)

    def resumen(self, mostrar=True):
        lineas = [f"(n,m) = ({self.n},{self.m}), k = {self.k}, "
                  f"b(n,m) = {self.maximo_puentes}, candidatos = {self.candidatos}",
                  f"  F-tupla maxima: {list(self.tupla_maxima)}",
                  f"  decide en F_{self.indice_discriminante}; "
                  f"clases ganadoras no isomorfas: {len(self.ganadores)}"]
        for H in self.ganadores:
            d = H.descripcion()
            lineas.append(f"  - aristas {sorted(H.grafo().edges(labels=False))}, "
                          f"terminales {list(H.terminales())}")
            lineas.append(f"    puentes {d['puentes']}, nucleos {d['nucleos']}")
        texto = "\n".join(lineas)
        if mostrar:
            print(texto)
        return texto


def ganadores(n, m, k, verbose=False):
    r"""
    Grafos k-terminales en T_{n,m}^{(k)} con F-tupla lexicograficamente maxima
    (los localmente mas separables cerca de p = 1).
    """
    mejor, segunda = None, None
    lista = []
    total = 0
    maximo_puentes = 0
    for G in grafos_conexos(n, m):
        maximo_puentes = max(maximo_puentes, len(list(G.bridges(labels=False))))
        if G.order() < k:
            continue
        for H in candidatos(G, k):
            total += 1
            t = H.F_tupla()
            if mejor is None or t > mejor:
                if mejor is not None:
                    segunda = mejor
                mejor, lista = t, [H]
            elif t == mejor:
                lista.append(H)
            elif segunda is None or t > segunda:
                segunda = t
    if mejor is None:
        raise ValueError(f"T_{{{n},{m}}} vacio para k = {k}")
    r = ResultadoFiltro(n, m, k, mejor, lista, segunda, maximo_puentes, total)
    if verbose:
        r.resumen()
    return r


def mejor_con_restriccion(n, m, k, condicion):
    r"""
    Mejor F-tupla entre los grafos k-terminales de T_{n,m} que cumplen
    condicion(H) (por ejemplo, tener al menos un puente).
    Devuelve (tupla, lista de GrafoKTerminales) o (None, []).
    """
    mejor, lista = None, []
    for G in grafos_conexos(n, m):
        if G.order() < k:
            continue
        for H in candidatos(G, k):
            if not condicion(H):
                continue
            t = H.F_tupla()
            if mejor is None or t > mejor:
                mejor, lista = t, [H]
            elif t == mejor:
                lista.append(H)
    return mejor, lista


# ----------------------------------------------------------------------
# Globos (Romero 2025, Definiciones 8 y 9)
# ----------------------------------------------------------------------

def en_I0(n, m):
    return n >= 4 and binomial(n - 1, 2) + 2 <= m <= binomial(n, 2)


def en_I(n, m):
    return n >= 4 and n <= m <= binomial(n, 2)


def grafo_globo(n, m):
    r"""Globo B_{n,m} (Definicion 8), definido recursivamente para (n, m) en I."""
    if not en_I(n, m):
        raise ValueError(f"({n},{m}) no esta en I")
    if en_I0(n, m):
        B = graphs.CompleteGraph(n - 1)
        B.add_vertex(n - 1)
        for u in range(m - binomial(n - 1, 2)):
            B.add_edge(n - 1, u)
        return B
    if (n, m) == (4, 4):
        B = graphs.CompleteGraph(3)
        B.add_edge(0, 3)
        return B
    B = grafo_globo(n - 1, m - 1)
    grados = B.degree(labels=True)
    v = min(B.vertices(sort=True), key=lambda x: (grados[x], x))
    B.add_edge(v, n - 1)
    return B


def globo_dos_terminales(n, m):
    r"""G_{n,m}: el globo con dos terminales a distancia igual a su diametro (Definicion 9)."""
    B = grafo_globo(n, m)
    D = B.distance_all_pairs()
    diametro = B.diameter()
    s, t = next((u, v) for u, v in combinations(B.vertices(sort=True), 2)
                if D[u][v] == diametro)
    return GrafoKTerminales(B, [s, t])


def verificar_teorema10(n_max=8, n_min=4, verbose=True):
    r"""
    Para cada (n, m) en I con n_min <= n <= n_max, comprueba que la F-tupla
    del globo G_{n,m} es la maxima de T_{n,m} (k = 2). Devuelve las fallas.
    """
    fallas = []
    for n in range(n_min, n_max + 1):
        for m in range(n, binomial(n, 2) + 1):
            r = ganadores(n, m, 2)
            G = globo_dos_terminales(n, m)
            ok = G.F_tupla() == r.tupla_maxima
            if not ok:
                fallas.append((n, m))
            if verbose:
                print(f"({n},{m}): {'ok' if ok else 'FALLA'}  "
                      f"decide en F_{r.indice_discriminante}, "
                      f"clases ganadoras: {len(r.ganadores)}, "
                      f"candidatos: {r.candidatos}")
    if verbose:
        print("Teorema 10 reproducido" if not fallas else f"Fallas: {fallas}")
    return fallas


# ----------------------------------------------------------------------
# Datos de k = 3 (exploracion inicial, n = 6, 7)
# ----------------------------------------------------------------------

# Prefijo (F_1, F_2, ...) de la F-tupla maxima hasta el indice donde quedo
# una sola clase, segun la primera corrida exhaustiva.
PREFIJOS_K3 = {
    (6, 6): [0, 3, 16], (6, 7): [0, 1, 9], (6, 8): [0, 1, 6, 19, 35, 0],
    (6, 9): [0, 0, 1, 7, 26], (6, 10): [0, 0, 0, 2, 12, 39],
    (6, 11): [0, 0, 0, 0, 2, 13], (6, 12): [0, 0, 0, 0, 0, 3],
    (7, 7): [0, 5, 24], (7, 8): [0, 3, 19, 45], (7, 9): [0, 3, 18, 51],
    (7, 10): [0, 1, 10, 41, 97], (7, 11): [0, 1, 9, 40],
    (7, 12): [0, 1, 10, 45, 124, 229], (7, 13): [0, 0, 1, 10, 46, 132],
    (7, 14): [0, 0, 0, 1, 11, 54], (7, 15): [0, 0, 0, 0, 2, 20],
    (7, 16): [0, 0, 0, 0, 0, 2],
}


def tabla_k3(ns=(6, 7), verbose=True):
    r"""
    Ganadores para k = 3 en cada (n, m) con n en ns y m >= n. Compara con
    PREFIJOS_K3 y marca si el esqueleto del ganador es casi-completo
    (la Proposicion 17 de Wurth). Devuelve la lista de ResultadoFiltro.
    """
    resultados = []
    for n in ns:
        for m in range(n, binomial(n, 2) + 1):
            r = ganadores(n, m, 3)
            resultados.append(r)
            if not verbose:
                continue
            prefijo = PREFIJOS_K3.get((n, m))
            control = ""
            if prefijo is not None:
                coincide = list(r.tupla_maxima[:len(prefijo)]) == prefijo
                control = "  prefijo ok" if coincide else "  PREFIJO DISTINTO"
            wurth = all(H.esqueleto_casi_completo() for H in r.ganadores)
            print(f"({n},{m}) b = {r.maximo_puentes}: decide en "
                  f"F_{r.indice_discriminante}, clases {len(r.ganadores)}, "
                  f"esqueleto casi-completo: {wurth}{control}")
            for H in r.ganadores:
                print("    ", H.descripcion()["nucleos"])
    return resultados


# ----------------------------------------------------------------------
# Arboles: F_{k-1} como menor del Laplaciano
# ----------------------------------------------------------------------

def verificar_menor_arbol(n_max=10, ks=(3, 4, 5), verbose=True):
    r"""
    Comprueba F_{k-1}(T) = det L(T)_[R] para todo arbol T con n <= n_max
    vertices y todo conjunto R de k terminales. Devuelve (chequeos, fallas).
    """
    chequeos, fallas = 0, []
    for n in range(2, n_max + 1):
        for T in graphs.trees(n):
            tabla = TablaConexos(T)
            for k in ks:
                if k > n:
                    continue
                for R in combinations(T.vertices(sort=True), k):
                    H = GrafoKTerminales(T, R, tabla=tabla)
                    chequeos += 1
                    if H.F(k - 1) != H.menor_laplaciano_terminales():
                        fallas.append((T.graph6_string(), R))
    if verbose:
        print(f"{chequeos} pares (arbol, terminales); fallas: {len(fallas)}")
    return chequeos, fallas


def tau(k, L):
    r"""
    Maximo de F_{k-1} sobre arboles con L aristas cuyas hojas son exactamente
    los k terminales, y los arboles que lo alcanzan (via el menor del Laplaciano).
    """
    mejor, lista = -1, []
    for T in graphs.trees(L + 1):
        hojas = [v for v in T if T.degree(v) == 1]
        if len(hojas) != k:
            continue
        valor = GrafoKTerminales(T, hojas).menor_laplaciano_terminales()
        if valor > mejor:
            mejor, lista = valor, [T]
        elif valor == mejor:
            lista.append(T)
    return mejor, lista

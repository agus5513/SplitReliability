r"""
GrafoKTerminales.sage

Grafos simples conexos con k terminales y su separabilidad.

Un subgrafo generador de G es k-separable si tiene exactamente k componentes
conexas, cada una con exactamente un terminal. Para i = 0, ..., m,

    F_i(G) = # subgrafos k-separables que se obtienen borrando i aristas,

    S_G(p) = sum_i F_i(G) (1 - p)^i p^(m - i).

G es localmente mas separable cerca de p = 1 si y solo si su F-tupla
(F_1, ..., F_m) es lexicograficamente maxima en su clase T_{n,m}^{(k)}.

Dos metodos independientes para la F-tupla (aritmetica entera exacta):

  - "particiones" (por defecto): suma sobre particiones de V en k bloques,
    un terminal por bloque, del producto de los polinomios de subgrafos
    conexos generadores de cada bloque. Costo ~ 3^n + k^(n-k) por grafo.
  - "aristas": fuerza bruta sobre subconjuntos de aristas borradas, con
    union-find. Sirve como control cruzado.

Uso:
    load("DataTypes/GrafoKTerminales.sage")
    G = GrafoKTerminales(graphs.CompleteGraph(4), [0, 1, 2])
    G.F_tupla()
    G.polinomio_separabilidad()
"""

from itertools import combinations, product


def _sumar_pol(a, b):
    r"""Suma de polinomios dados como listas de coeficientes."""
    if len(a) < len(b):
        a, b = b, a
    r = list(a)
    for i, c in enumerate(b):
        r[i] += c
    return r


def _multiplicar_pol(a, b):
    r"""Producto de polinomios dados como listas de coeficientes."""
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    r[i + j] += x * y
    return r


def _restar_pol(a, b):
    return _sumar_pol(a, [-c for c in b])


def es_casi_completo(H):
    r"""True si H tiene al menos binomial(n-1, 2) + 2 aristas (Wurth, Def. en Sec. 1)."""
    n = H.order()
    return H.size() >= binomial(n - 1, 2) + 2


def maximo_puentes_formula(n, m):
    r"""
    b(n, m) por la formula corregida, valida para n <= m <= binomial(n, 2):
    b(n, m) = n - 1 - floor((3 + sqrt(8(m - n) + 1)) / 2).
    Para trabajo numerico preferir el conteo directo (ver FiltroLexicografico).
    """
    if m == n - 1:
        return n - 1
    return n - 1 - (3 + Integer(8 * (m - n) + 1).isqrt()) // 2


class TablaConexos:
    r"""
    Polinomios de subgrafos conexos generadores de cada subgrafo inducido G[S].

    Para S subconjunto de V (como mascara de bits), C(S) es la lista cuyo
    coeficiente j cuenta los subgrafos generadores conexos de G[S] con
    exactamente j aristas. Se calcula con la recurrencia clasica

        C(S) = A(S) - sum_{v0 in T, T subconjunto propio de S} C(T) A(S \ T),

    donde A(S) = (1 + y)^{e(S)} cuenta todos los subgrafos de G[S] y v0 es el
    menor vertice de S. Se memoiza por mascara: una misma tabla se comparte
    entre todos los conjuntos de terminales de un mismo grafo.
    """

    def __init__(self, G):
        self.vertices = list(G.vertices(sort=True))
        self.indice = {v: i for i, v in enumerate(self.vertices)}
        n = len(self.vertices)
        self.adyacencia = [0] * n
        for u, v in G.edges(labels=False, sort=False):
            iu, iv = self.indice[u], self.indice[v]
            self.adyacencia[iu] |= 1 << iv
            self.adyacencia[iv] |= 1 << iu
        self._binomiales = [[binomial(e, j) for j in range(e + 1)]
                            for e in range(G.size() + 1)]
        self._C = {}
        self._aristas = {}

    def aristas(self, S):
        r"""Cantidad de aristas de G[S] (memoizada por mascara)."""
        if S not in self._aristas:
            self._aristas[S] = self._contar_aristas(S)
        return self._aristas[S]

    def _contar_aristas(self, S):
        total = 0
        T = S
        while T:
            v = (T & -T).bit_length() - 1
            total += bin(self.adyacencia[v] & S).count("1")
            T &= T - 1
        return total // 2

    def todos(self, S):
        r"""A(S) = (1 + y)^{e(S)}."""
        return self._binomiales[self.aristas(S)]

    def C(self, S):
        r"""Polinomio (por aristas conservadas) de subgrafos conexos generadores de G[S]."""
        if S in self._C:
            return self._C[S]
        if S & (S - 1) == 0:
            resultado = [1]
        else:
            v0 = S & -S
            resto = S - v0  # v0 esta en S (no usar ^ ni ~: Sage los reinterpreta)
            resultado = self.todos(S)
            # T recorre los subconjuntos propios de S que contienen a v0.
            U = (resto - 1) & resto
            while True:
                T = v0 | U
                if T != S:
                    resultado = _restar_pol(
                        resultado, _multiplicar_pol(self.C(T), self.todos(S - T)))
                if U == 0:
                    break
                U = (U - 1) & resto
        while len(resultado) > 1 and resultado[-1] == 0:
            resultado.pop()
        self._C[S] = resultado
        return resultado


class GrafoKTerminales:
    r"""
    Grafo simple conexo con k terminales distinguidos.

    INPUT:
        - G: un Graph de Sage (o algo que Graph() acepte), simple y conexo.
        - terminales: lista de k vertices distintos de G.
        - tabla: (opcional) TablaConexos de G, para compartirla entre
          distintos conjuntos de terminales del mismo grafo.
    """

    def __init__(self, G, terminales, tabla=None):
        G = Graph(G, immutable=True)
        if G.has_multiple_edges() or G.has_loops():
            raise ValueError("el grafo tiene que ser simple")
        if not G.is_connected():
            raise ValueError("el grafo tiene que ser conexo")
        terminales = tuple(terminales)
        if len(set(terminales)) != len(terminales):
            raise ValueError("terminales repetidos")
        if any(t not in G for t in terminales):
            raise ValueError("hay terminales que no son vertices de G")
        self._G = G
        self._terminales = terminales
        self._tabla = tabla if tabla is not None else TablaConexos(G)
        self._F = None

    # ------------------------------------------------------------------
    # Datos basicos
    # ------------------------------------------------------------------

    def grafo(self):
        return self._G

    def terminales(self):
        return self._terminales

    def n(self):
        return self._G.order()

    def m(self):
        return self._G.size()

    def k(self):
        return len(self._terminales)

    def __repr__(self):
        return (f"Grafo {self.k()}-terminal con n = {self.n()}, m = {self.m()}, "
                f"terminales {list(self._terminales)}")

    # ------------------------------------------------------------------
    # Estructura: puentes, esqueleto, nucleos
    # ------------------------------------------------------------------

    def puentes(self):
        return list(self._G.bridges(labels=False))

    def numero_puentes(self):
        return len(self.puentes())

    def esqueleto(self):
        r"""Grafo que resulta de contraer todos los puentes."""
        H = self._G.copy(immutable=False)
        H.contract_edges(self.puentes())
        return H

    def nucleos(self):
        r"""Componentes 2-arista-conexas con mas de un vertice, como subgrafos inducidos."""
        H = self._G.copy(immutable=False)
        H.delete_edges(self.puentes())
        return [self._G.subgraph(c) for c in H.connected_components()
                if len(c) > 1]

    def esqueleto_casi_completo(self):
        return es_casi_completo(self.esqueleto())

    def descripcion(self):
        r"""
        Resumen estructural para leer tablas: puentes, nucleo, grados del nucleo,
        terminales dentro del nucleo y vertices de enganche de los puentes.
        """
        puentes = self.puentes()
        nucleos = self.nucleos()
        resumen = {"n": self.n(), "m": self.m(), "puentes": len(puentes),
                   "nucleos": []}
        for N in nucleos:
            grados = N.degree(labels=True)
            enganches = sorted({x for e in puentes for x in e if x in N})
            resumen["nucleos"].append({
                "n'": N.order(),
                "m'": N.size(),
                "casi_completo": es_casi_completo(N),
                "grados": sorted(grados.values()),
                "terminales_dentro": [(t, grados[t]) for t in self._terminales
                                      if t in N],
                "enganches": [(x, grados[x]) for x in enganches],
            })
        return resumen

    # ------------------------------------------------------------------
    # Separabilidad
    # ------------------------------------------------------------------

    def F_polinomio(self, metodo="particiones"):
        r"""
        Lista [F_0, F_1, ..., F_m]. metodo: "particiones" (rapido) o "aristas"
        (fuerza bruta, solo para controles en grafos chicos).
        """
        if metodo == "aristas":
            return self._F_por_aristas()
        if metodo != "particiones":
            raise ValueError("metodo desconocido: " + str(metodo))
        if self._F is None:
            self._F = self._F_por_particiones()
        return list(self._F)

    def F(self, i):
        return self.F_polinomio()[i]

    def F_tupla(self, metodo="particiones"):
        r"""(F_1, ..., F_m), la tupla que se compara lexicograficamente."""
        return tuple(self.F_polinomio(metodo)[1:])

    def primer_coeficiente_no_nulo(self):
        r"""Menor i con F_i > 0 (es al menos k - 1)."""
        F = self.F_polinomio()
        return next(i for i, c in enumerate(F) if c)

    def polinomio_separabilidad(self, anillo=None):
        r"""S_G(p) como elemento de QQ[p]."""
        if anillo is None:
            anillo = PolynomialRing(QQ, "p")
        p = anillo.gen()
        m = self.m()
        return sum(c * (1 - p)**i * p**(m - i)
                   for i, c in enumerate(self.F_polinomio()))

    def _F_por_particiones(self):
        tabla = self._tabla
        n, m = self.n(), self.m()
        ts = [tabla.indice[t] for t in self._terminales]
        otros = [i for i in range(n) if i not in ts]
        K = [0]
        for asignacion in product(range(len(ts)), repeat=len(otros)):
            bloques = [1 << t for t in ts]
            for v, a in zip(otros, asignacion):
                bloques[a] |= 1 << v
            termino = [1]
            for B in bloques:
                termino = _multiplicar_pol(termino, tabla.C(B))
                if not any(termino):
                    break
            K = _sumar_pol(K, termino)
        # K[j] = # subgrafos k-separables con j aristas conservadas.
        F = [0] * (m + 1)
        for j, c in enumerate(K):
            if c:
                F[m - j] += c
        return F

    def _F_por_aristas(self):
        aristas = list(self._G.edges(labels=False, sort=True))
        m, k = len(aristas), self.k()
        indice = {v: i for i, v in enumerate(self._G.vertices(sort=True))}
        n = len(indice)
        ts = [indice[t] for t in self._terminales]
        pares = [(indice[u], indice[v]) for u, v in aristas]
        F = [0] * (m + 1)
        for i in range(m + 1):
            for borradas in combinations(range(m), i):
                fuera = set(borradas)
                padre = list(range(n))

                def raiz(x):
                    while padre[x] != x:
                        padre[x] = padre[padre[x]]
                        x = padre[x]
                    return x

                for j, (u, v) in enumerate(pares):
                    if j not in fuera:
                        padre[raiz(u)] = raiz(v)
                raices = {raiz(v) for v in range(n)}
                if len(raices) == k and len({raiz(t) for t in ts}) == k:
                    F[i] += 1
        return F

    # ------------------------------------------------------------------
    # Arboles: F_{k-1} como menor del Laplaciano
    # ------------------------------------------------------------------

    def menor_laplaciano_terminales(self):
        r"""det de L(G) sin las filas y columnas de los terminales (bosques enraizados)."""
        V = list(self._G.vertices(sort=True))
        L = self._G.laplacian_matrix(vertices=V)
        quitar = [V.index(t) for t in self._terminales]
        M = L.delete_rows(quitar).delete_columns(quitar)
        return M.det() if M.nrows() else 1

    # ------------------------------------------------------------------
    # Isomorfismo con terminales (los terminales pueden permutarse entre si)
    # ------------------------------------------------------------------

    def forma_canonica(self):
        r"""Clave hashable: dos grafos k-terminales son isomorfos sii coinciden."""
        terminales = list(self._terminales)
        resto = [v for v in self._G.vertices(sort=True) if v not in terminales]
        particion = [terminales, resto] if resto else [terminales]
        canonico, certificado = self._G.canonical_label(partition=particion,
                                                       certificate=True)
        return (canonico.graph6_string(),
                tuple(sorted(certificado[t] for t in terminales)))

    def es_isomorfo(self, otro):
        return self.forma_canonica() == otro.forma_canonica()

    def es_equivalente(self, otro):
        r"""Misma F-tupla, es decir, mismo polinomio de separabilidad (Remark 5)."""
        return self.F_tupla() == otro.F_tupla()

    def dibujar(self, **opciones):
        r"""Plot con los terminales resaltados y los puentes en rojo."""
        colores = {"lightblue": list(self._terminales)}
        return self._G.plot(vertex_colors=colores,
                            edge_colors={"red": self.puentes()}, **opciones)

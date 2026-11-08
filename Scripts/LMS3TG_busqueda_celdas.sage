# %% [markdown]
# # Búsqueda exhaustiva: los LMS3TG son el grafo modelo $G^*_{n,m}$
#
# Para cada clase $\mathcal{T}_{n,m}$ con $b(n,m)\ge 2$ este notebook:
#
# 1. enumera todos los grafos 3-terminales no isomorfos (nauty, más una órbita de terminales por grafo);
# 2. calcula su tupla $F=(F_1,\dots,F_m)$ con aritmética entera exacta;
# 3. se queda con los de tupla lexicográficamente máxima (los localmente más separables cerca de $p=1$);
# 4. construye $G^*_{n,m}$ según la Definición del grafo modelo y comprueba que tiene esa tupla, que los ganadores son exactamente su clase y que el índice de decisión es el predicho.
#
# Es autocontenido: no usa los módulos `.sage` anteriores. Necesita SageMath con `nauty` (incluido en Sage).

# %% [markdown]
# ## 1. Separabilidad: la tupla $F$
#
# $F_i(G)$ se calcula sumando, sobre las particiones de $V$ en tres bloques con un terminal cada uno, el producto de los polinomios de subgrafos conexos generadores de cada bloque. Esos polinomios salen de la recurrencia estándar sobre subconjuntos (máscaras de bits). Nota: en Sage `^` es potencia, por eso las restas de máscaras se escriben `S - T`.

# %% [celda 1]
from itertools import combinations, product


def _suma(a, b):
    if len(a) < len(b):
        a, b = b, a
    r = list(a)
    for i, c in enumerate(b):
        r[i] += c
    return r


def _producto(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    r[i + j] += x * y
    return r


class TablaConexos:
    # C(S): coeficiente j = # subgrafos generadores conexos de G[S] con j aristas
    def __init__(self, G):
        self.vertices = list(G.vertices(sort=True))
        self.indice = {v: i for i, v in enumerate(self.vertices)}
        self.ady = [0] * len(self.vertices)
        for u, v in G.edges(labels=False, sort=False):
            iu, iv = self.indice[u], self.indice[v]
            self.ady[iu] |= 1 << iv
            self.ady[iv] |= 1 << iu
        self.binom = [[binomial(e, j) for j in range(e + 1)] for e in range(G.size() + 1)]
        self._C, self._e = {}, {}

    def aristas(self, S):
        if S not in self._e:
            total, T = 0, S
            while T:
                v = (T & -T).bit_length() - 1
                total += bin(self.ady[v] & S).count("1")
                T &= T - 1
            self._e[S] = total // 2
        return self._e[S]

    def todos(self, S):
        return self.binom[self.aristas(S)]

    def C(self, S):
        if S in self._C:
            return self._C[S]
        if S & (S - 1) == 0:
            res = [1]
        else:
            v0 = S & -S
            resto = S - v0
            res = list(self.todos(S))
            U = (resto - 1) & resto
            while True:
                T = v0 | U
                if T != S:
                    res = _suma(res, [-c for c in _producto(self.C(T), self.todos(S - T))])
                if U == 0:
                    break
                U = (U - 1) & resto
        while len(res) > 1 and res[-1] == 0:
            res.pop()
        self._C[S] = res
        return res


def F_tupla(G, terminales, tabla=None):
    # (F_1, ..., F_m) de G con los terminales dados
    tabla = tabla or TablaConexos(G)
    n, m = G.order(), G.size()
    ts = [tabla.indice[t] for t in terminales]
    otros = [i for i in range(n) if i not in ts]
    K = [0]
    for asign in product(range(3), repeat=len(otros)):
        bloques = [1 << t for t in ts]
        for v, a in zip(otros, asign):
            bloques[a] |= 1 << v
        term = [1]
        for B in bloques:
            term = _producto(term, tabla.C(B))
            if not any(term):
                break
        K = _suma(K, term)
    F = [0] * (m + 1)
    for j, c in enumerate(K):          # j aristas conservadas = m - j borradas
        F[m - j] += c
    return tuple(F[1:])


def forma_canonica(G, terminales):
    # clave de isomorfismo de grafos 3-terminales (los terminales pueden permutarse)
    ts = list(terminales)
    resto = [v for v in G.vertices(sort=True) if v not in ts]
    can, cert = G.canonical_label(partition=[ts, resto] if resto else [ts], certificate=True)
    return (can.graph6_string(), tuple(sorted(cert[t] for t in ts)))


# control rapido: K_{1,3} con terminales en las hojas tiene F = (0, 3, 0)
print(F_tupla(Graph([(0, 1), (0, 2), (0, 3)]), [1, 2, 3]))

# %% [markdown]
# ## 2. El grafo modelo $G^*_{n,m}$
#
# Implementa la Definición: $b=b(n,m)$, $n'=n-b$, $m'=m-b$, $\lambda'=m'-\binom{n'-1}{2}$, $B=B_{n',m'}$ y la terna justa $\ell(b)$. El globo se construye con $K_{n'-1}$ en los vértices $0,\dots,n'-2$ y $z=n'-1$ adyacente a $0,\dots,\lambda'-1$. Así, en el caso $\lambda'=n'-2$ los vértices de grado mínimo son $p=z$ y $q=n'-2$, y en el caso $\lambda'\le n'-3$ los no vecinos de $z$ son $W=\{\lambda',\dots,n'-2\}$.

# %% [celda 2]
def b_nm(n, m):
    # b(n,m) = max{b : m - b <= binom(n-b, 2)}, para n <= m <= binom(n,2)
    return max(b for b in range(n) if m - b <= binomial(n - b, 2))


def terna_justa(b):
    q, r = divmod(b, 3)
    return tuple(sorted([q + (1 if i < r else 0) for i in range(3)], reverse=True))


def globo(np_, mp_):
    lam = mp_ - binomial(np_ - 1, 2)
    B = graphs.CompleteGraph(np_ - 1)
    B.add_vertex(np_ - 1)
    for u in range(lam):
        B.add_edge(np_ - 1, u)
    return B, np_ - 1, lam


def theta(H, patas):
    # Theta(H; (u1,l1), (u2,l2), (u3,l3)): devuelve (grafo, terminales)
    G = H.copy(immutable=False)
    sig = max(G.vertices(sort=True)) + 1
    ts = []
    for u, l in patas:
        actual = u
        for _ in range(l):
            G.add_edge(actual, sig)
            actual, sig = sig, sig + 1
        ts.append(actual)
    return G, ts


def caso(n, m):
    b = b_nm(n, m)
    np_, mp_ = n - b, m - b
    lam = mp_ - binomial(np_ - 1, 2)
    l1, l2, l3 = terna_justa(b)
    if lam == np_ - 1:
        return "(a)"
    if lam == np_ - 2:
        return "(b)" if l2 == l3 else "(b')"
    return "(c)"


def modelo(n, m):
    b = b_nm(n, m)
    np_, mp_ = n - b, m - b
    B, z, lam = globo(np_, mp_)
    l1, l2, l3 = terna_justa(b)
    c = caso(n, m)
    if c == "(a)":
        return theta(B, [(0, l1), (1, l2), (2, l3)])
    p, q = z, np_ - 2
    if c == "(b)":
        return theta(B, [(0, l1), (p, l2), (q, l3)])
    if c == "(b')":
        return theta(B, [(p, l1), (p, l2), (q, l3)])
    return theta(B, [(lam, l1), (lam + 1, l2), (z, l3)])   # w1, w2 en W


def variantes_bprima(n, m):
    # caso (b') con b >= 5: nucleo sobre la pata corta (en su extremo o en el interior)
    b = b_nm(n, m)
    if caso(n, m) != "(b')" or b < 3:
        return []
    np_, mp_ = n - b, m - b
    B, z, lam = globo(np_, mp_)
    l1, l2, l3 = terna_justa(b)
    p, q = z, np_ - 2
    salida = []
    for d in range(1, l3 + 1):          # distancia del nucleo a la mediana
        G = B.copy(immutable=False)
        sig = max(G.vertices(sort=True)) + 1
        actual = q
        for _ in range(d):
            G.add_edge(actual, sig); actual, sig = sig, sig + 1
        mediana, ts = actual, []
        for l in (l1, l2):
            actual = mediana
            for _ in range(l):
                G.add_edge(actual, sig); actual, sig = sig, sig + 1
            ts.append(actual)
        actual = p
        for _ in range(l3 - d):
            G.add_edge(actual, sig); actual, sig = sig, sig + 1
        ts.append(actual)
        salida.append((G, ts))
    return salida


def indice_predicho(n, m):
    # Tabla del Teorema A
    b = b_nm(n, m)
    np_ = n - b
    lam = m - b - binomial(np_ - 1, 2)
    c = caso(n, m)
    if c == "(a)":
        return np_ if b >= 3 else 2 * np_ - 3
    if c == "(b)":
        return np_
    if c == "(b')":
        return np_ - 1
    return np_ - 1 if b >= 3 else lam + np_ - 2


G, ts = modelo(7, 10)
print("G*_{7,10}:", sorted(G.edges(labels=False)), "terminales", ts)
print("F =", F_tupla(G, ts))

# %% [markdown]
# ## 3. Búsqueda exhaustiva en $\mathcal{T}_{n,m}$

# %% [celda 3]
def ganadores(n, m):
    # devuelve (tupla maxima, lista de ganadores no isomorfos, segunda mejor tupla, #candidatos)
    mejor, segunda, lista, total = None, None, [], 0
    for G in graphs.nauty_geng(f"{n} {m}:{m} -c"):
        tabla = TablaConexos(G)
        vistos = set()
        for ts in combinations(G.vertices(sort=True), 3):
            clave = forma_canonica(G, ts)
            if clave in vistos:
                continue
            vistos.add(clave)
            total += 1
            t = F_tupla(G, ts, tabla)
            if mejor is None or t > mejor:
                if mejor is not None:
                    segunda = mejor
                mejor, lista = t, [(G, ts)]
            elif t == mejor:
                lista.append((G, ts))
            elif segunda is None or t > segunda:
                segunda = t
    return mejor, lista, segunda, total


def primer_indice_distinto(a, b):
    return next((i + 1 for i, (x, y) in enumerate(zip(a, b)) if x != y), None)

# %% [markdown]
# ## 4. Verificación del Teorema A
#
# Para cada $(n,m)$ con $n\le m\le\binom{n-2}{2}+2$ (es decir, $b(n,m)\ge2$) se comprueba:
#
# - la tupla máxima es la de $G^*_{n,m}$;
# - cada ganador es isomorfo a $G^*_{n,m}$ o, en el caso (b'), a una de sus variantes;
# - el índice de decisión observado (primer índice donde la mejor tupla le gana a la segunda) es el de la tabla del teorema.
#
# Tiempos aproximados: $n\le7$ en menos de un minuto; $n=8$ del orden de media hora.

# %% [celda 4]
N_MAX = 7          # subir a 8 para la verificacion completa de n = 8

filas = []
for n in range(5, N_MAX + 1):
    for m in range(n, binomial(n - 2, 2) + 3):
        if b_nm(n, m) < 2:
            continue
        mejor, lista, segunda, total = ganadores(n, m)
        G, ts = modelo(n, m)
        modelos = [(G, ts)] + variantes_bprima(n, m)
        claves_modelo = {forma_canonica(H, t) for H, t in modelos}
        tupla_ok = all(F_tupla(H, t) == mejor for H, t in modelos)
        clases_ok = {forma_canonica(H, t) for H, t in lista} == claves_modelo
        d_obs = primer_indice_distinto(mejor, segunda)
        d_pred = indice_predicho(n, m)
        b = b_nm(n, m)
        np_ = n - b
        lam = m - b - binomial(np_ - 1, 2)
        filas.append((n, m, b, np_, lam, caso(n, m), terna_justa(b), total, len(lista),
                      tupla_ok, clases_ok, d_pred, d_obs))
        print(f"({n},{m}) b={b} n'={np_} lam'={lam} caso {caso(n, m):4} l={terna_justa(b)} "
              f"candidatos={total:6}  ganadores={len(lista)}  "
              f"tupla={'OK' if tupla_ok else 'FALLA'}  clases={'OK' if clases_ok else 'FALLA'}  "
              f"indice pred/obs = {d_pred}/{d_obs}", flush=True)

todo_ok = all(f[9] and f[10] and f[11] == f[12] for f in filas)
print()
print("Todas las clases coinciden con G*_{n,m} y con el indice predicho:" if todo_ok else "HAY DISCREPANCIAS", todo_ok)

# %% [markdown]
# ## 5. Un caso en detalle: dónde pierde el segundo mejor

# %% [celda 5]
n, m = 7, 10
mejor, lista, segunda, total = ganadores(n, m)
d = primer_indice_distinto(mejor, segunda)
print("tupla de G*  :", mejor)
print("segunda mejor:", segunda)
print("deciden en F_%d: %d contra %d" % (d, mejor[d - 1], segunda[d - 1]))
G, ts = modelo(n, m)
G.plot(vertex_colors={"lightblue": ts}, edge_colors={"red": list(G.bridges(labels=False))}, figsize=4)

# %% [markdown]
# ## 6. El caso (b'): una clase con grafos no isomorfos
#
# En $(9,10)$ se tiene $b=5$, $\ell=(2,2,1)$ y núcleo $K_4-pq$. El teorema predice dos ganadores no isomorfos con la misma tupla: $G^*$ (dos patas largas en $p$, la corta en $q$) y el grafo con el núcleo en el extremo de la pata corta. Tarda unos minutos.

# %% [celda 6]
mejor, lista, segunda, total = ganadores(9, 10)
print("candidatos:", total, "| ganadores no isomorfos:", len(lista), "| tupla:", mejor)
modelos = [modelo(9, 10)] + variantes_bprima(9, 10)
print("son exactamente G* y su variante:",
      {forma_canonica(H, t) for H, t in lista} == {forma_canonica(H, t) for H, t in modelos})
graphics_array([H.plot(vertex_colors={"lightblue": t},
                       edge_colors={"red": list(H.bridges(labels=False))}) for H, t in lista]).show(figsize=[8, 3])

# %% [markdown]
# ## 7. Fuera del teorema: $b(n,m)=1$
#
# Con $m=\binom{n-1}{2}+1$ la plantilla $G^*_{n,m}$ (con terna $(1,0,0)$) no es el ganador: el LMS3TG no tiene puentes.

# %% [celda 7]
for n in (6, 7):
    m = binomial(n - 1, 2) + 1
    mejor, lista, segunda, total = ganadores(n, m)
    G, ts = modelo(n, m)
    H, th = lista[0]
    print(f"({n},{m}) b={b_nm(n, m)}: puentes del ganador = {len(list(H.bridges()))}, "
          f"F(ganador) = {mejor[:2*n-5]}..., F(G*) = {F_tupla(G, ts)[:2*n-5]}...")

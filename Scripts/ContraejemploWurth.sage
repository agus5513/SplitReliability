r"""
ContraejemploWurth.sage

Certificado del contraejemplo a la Proposicion 17 de Wurth y Romero (2026)
("el esqueleto de todo LMS3TG es casi-completo") cuando b(n, m) = 1.

Para (n, m) = (6, 11) y (7, 16), con k = 3:

  1. El grafo G de CONTRAEJEMPLOS no tiene puentes y G no es casi-completo.
  2. Su F-tupla coincide por dos metodos independientes (particiones y aristas).
  3. Busqueda exhaustiva: G tiene la F-tupla maxima de T_{n,m}, y le gana
     al mejor grafo con al menos un puente.
  4. Prueba exacta de que S_G(p) > S_H(p) cerca de p = 1:
        S_G - S_H = (1 - p)^j q(p),  con q(1) = F_j(G) - F_j(H) > 0.

Requiere DataTypes/GrafoKTerminales.sage y Scripts/FiltroLexicografico.sage.

Uso:
    certificado = certificar(6, 11)
    print(informe_latex(certificado))
"""

try:
    GrafoKTerminales, ganadores
except NameError:
    raise ImportError("Carga primero DataTypes/GrafoKTerminales.sage y "
                      "Scripts/FiltroLexicografico.sage")


CONTRAEJEMPLOS = {
    (6, 11): ([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (1, 4), (1, 5),
               (2, 3), (2, 4), (2, 5), (3, 5)], [0, 4, 5]),
    # Grafo umbral: K4 sobre {3,4,5,6}; 0 y 1 adyacentes a todo el K4; 2 a {5,6}.
    (7, 16): ([(0, 3), (0, 4), (0, 5), (0, 6), (1, 3), (1, 4), (1, 5), (1, 6),
               (2, 5), (2, 6), (3, 4), (3, 5), (3, 6), (4, 5), (4, 6), (5, 6)],
              [0, 1, 2]),
}


def contraejemplo(n, m):
    aristas, terminales = CONTRAEJEMPLOS[(n, m)]
    return GrafoKTerminales(Graph(aristas), terminales)


def primer_indice_distinto(a, b):
    r"""Menor i (desde 1) con a_i != b_i para F-tuplas (F_1, ..., F_m)."""
    return next(i + 1 for i, (x, y) in enumerate(zip(a, b)) if x != y)


def comparacion_exacta(G, H):
    r"""
    Escribe S_G - S_H = (1 - p)^j q(p) con j el primer indice donde difieren
    las F-tuplas. Devuelve (j, q, q(1)); q(1) = F_j(G) - F_j(H).
    """
    R = PolynomialRing(QQ, "p")
    p = R.gen()
    D = G.polinomio_separabilidad(R) - H.polinomio_separabilidad(R)
    j = primer_indice_distinto(G.F_tupla(), H.F_tupla())
    q, resto = D.quo_rem((1 - p)**j)
    if resto != 0:
        raise ArithmeticError("(1-p)^j no divide a S_G - S_H")
    return j, q, q(1)


def certificar(n, m, verbose=True):
    r"""Corre los cuatro controles del encabezado y devuelve un diccionario."""
    G = contraejemplo(n, m)
    c = {"n": n, "m": m, "G": G}

    # 1. Estructura
    c["puentes_G"] = G.numero_puentes()
    c["G_casi_completo"] = G.esqueleto_casi_completo()
    c["umbral_casi_completo"] = binomial(n - 1, 2) + 2

    # 2. Dos metodos
    c["F_G"] = G.F_tupla()
    c["metodos_coinciden"] = G.F_tupla() == G.F_tupla("aristas")

    # 3. Busqueda exhaustiva
    r = ganadores(n, m, 3)
    c["b(n,m)"] = r.maximo_puentes
    c["G_es_ganador"] = any(H.es_isomorfo(G) for H in r.ganadores)
    c["clases_ganadoras"] = len(r.ganadores)
    t_puente, con_puente = mejor_con_restriccion(
        n, m, 3, lambda H: H.numero_puentes() >= 1)
    H = con_puente[0]
    c["H"] = H
    c["F_H"] = t_puente

    # 4. Comparacion exacta cerca de p = 1
    j, q, q1 = comparacion_exacta(G, H)
    c["j"], c["q"], c["q(1)"] = j, q, q1
    c["valido"] = (c["puentes_G"] == 0 and not c["G_casi_completo"]
                   and c["metodos_coinciden"] and c["G_es_ganador"]
                   and c["b(n,m)"] == 1 and q1 > 0)

    if verbose:
        print(f"(n,m) = ({n},{m}), k = 3, b(n,m) = {c['b(n,m)']}")
        print(f"  G: {c['puentes_G']} puentes, m = {m} < "
              f"{c['umbral_casi_completo']} = C({n - 1},2)+2, "
              f"casi-completo: {c['G_casi_completo']}")
        print(f"  F(G) por dos metodos coincide: {c['metodos_coinciden']}")
        print(f"  G es el ganador lexicografico: {c['G_es_ganador']} "
              f"(clases ganadoras: {c['clases_ganadoras']})")
        print(f"  F(G)              = {list(c['F_G'])}")
        print(f"  F(mejor con puente) = {list(c['F_H'])}")
        print(f"  S_G - S_H = (1-p)^{j} q(p),  q(1) = {q1}")
        print("  CONTRAEJEMPLO VALIDO" if c["valido"] else "  REVISAR")
    return c


def informe_latex(c):
    r"""Fragmento LaTeX con el contraejemplo, listo para pegar en una nota."""
    G, H = c["G"], c["H"]
    j = c["j"]
    hasta = max(i for i, x in enumerate(c["F_G"]) if x) + 1
    columnas = range(j - 2 if j > 2 else 1, hasta + 1)

    def fila(nombre, F):
        return nombre + " & " + " & ".join(str(F[i - 1]) for i in columnas) + r" \\"

    aristas = ", ".join(f"{u}{v}" for u, v in sorted(G.grafo().edges(labels=False)))
    return "\n".join([
        r"\begin{example}[$(n,m)=(%d,%d)$, $k=3$]" % (c["n"], c["m"]),
        r"Sea $G$ el grafo con aristas $%s$ y terminales $\{%s\}$." % (
            aristas, ",".join(map(str, G.terminales()))),
        r"Se tiene $b(%d,%d)=1$, $G$ no tiene puentes y $m=%d<\binom{%d}{2}+2=%d$,"
        % (c["n"], c["m"], c["m"], c["n"] - 1, c["umbral_casi_completo"]),
        r"por lo que su esqueleto no es casi-completo. Sin embargo $G$ es LMS3TG:",
        r"\begin{center}",
        r"\begin{tabular}{l" + "r" * len(columnas) + "}",
        " & " + " & ".join(f"$F_{{{i}}}$" for i in columnas) + r" \\ \hline",
        fila(r"$G$", c["F_G"]),
        fila(r"$H$ (mejor con un puente)", c["F_H"]),
        r"\end{tabular}",
        r"\end{center}",
        r"y $S_G(p)-S_H(p)=(1-p)^{%d}q(p)$ con $q(1)=%s>0$." % (j, c["q(1)"]),
        r"\end{example}",
    ])

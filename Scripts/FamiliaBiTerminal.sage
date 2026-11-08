load("DataTypes/Grafo2Terminales.sage")

from sage.all import Graph, graphs

def B(n, a, b, c):
    """
        Construye el grafo con 2 terminales B(n; a, b, c), con c <= b <= a.
        Consiste en un clique de q = n - 2 vértices (0..q-1) y dos
        terminales s = q, t = q + 1 por fuera del clique. El terminal s
        se conecta con a vértices del clique y t se conecta con c vecinos
        de s y con b - c vértices que no son vecinos de s. Así, s tiene
        a - c vecinos exclusivos, t tiene b - c y hay c vecinos comunes.
    """
    q = n - 2
    if not (0 <= c <= b <= a):
        raise ValueError(f"Se requiere c <= b <= a, se recibió a={a}, b={b}, c={c}")
    if a > q or a + b - c > q:
        raise ValueError(f"Parámetros incompatibles con n={n}: se requiere a + b - c <= n - 2")

    G = graphs.CompleteGraph(q)
    s, t = q, q + 1
    G.add_vertices([s, t])
    G.add_edges((s, v) for v in range(a))
    G.add_edges((t, v) for v in range(c))              # vecinos comunes
    G.add_edges((t, v) for v in range(a, a + b - c))   # vecinos exclusivos de t
    return Grafo2Terminales(G, s, t)

def Numero2BosquesBiTerminal(n, a, b, c):
    """
        Número t_2 de H = B(n; a, b, c) según la fórmula cerrada:
        q = n - 2, D = a + b, e = D - 2c, r = q - c - e,
        t_2(H) = (q+2)^c (q+1)^e q^r (1 - c/(q+2) - e/(q+1) - r/q).
    """
    q = n - 2
    D = a + b
    e = D - 2 * c
    r = q - c - e
    return (q + 2)**c * (q + 1)**e * q**r * (1 - QQ(c) / (q + 2) - QQ(e) / (q + 1) - QQ(r) / q)

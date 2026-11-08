from math import cos, sin, pi

NARANJA = "#d9822b"


def patas_y_etiquetas(n, m):
    # nucleo B, patas [(u_r, l_r)] y etiquetas de vertices especiales, segun la Definicion del grafo modelo
    b = b_nm(n, m)
    np_, mp_ = n - b, m - b
    B, z, lam = globo(np_, mp_)
    l1, l2, l3 = terna_justa(b)
    c = caso(n, m)
    if c == "(a)":
        return B, [(0, l1), (1, l2), (2, l3)], {}
    p, q = z, np_ - 2
    if c == "(b)":
        return B, [(0, l1), (p, l2), (q, l3)], {0: "y", p: "p", q: "q"}
    if c == "(b')":
        return B, [(p, l1), (p, l2), (q, l3)], {p: "p", q: "q"}
    return B, [(lam, l1), (lam + 1, l2), (z, l3)], {lam: "w_1", lam + 1: "w_2", z: "z"}


def dibujar_modelo(n, m, titulo=None, R=1.0, paso=0.75):
    B, patas, etiq = patas_y_etiquetas(n, m)
    G, ts = theta(B, patas)
    np_ = B.order()

    # vertices del nucleo sobre un circulo; los enganches donde los pone la figura
    enganches = []
    for u, l in patas:
        if u not in enganches:
            enganches.append(u)
    deseados = {1: [90], 2: [180, 0], 3: [90, 210, 330]}[len(enganches)]
    libres_ang = [(90 + 360 * i / np_) % 360 for i in range(np_)]
    ang = {}
    for u, a in zip(enganches, deseados):
        mejor = min(libres_ang, key=lambda x: min(abs(x - a), 360 - abs(x - a)))
        ang[u] = mejor
        libres_ang.remove(mejor)
    for v in B.vertices(sort=True):
        if v not in ang:
            ang[v] = libres_ang.pop(0)
    pos = {v: (R * cos(float(a) * pi / 180), R * sin(float(a) * pi / 180)) for v, a in ang.items()}

    # patas: caminos hacia afuera en la direccion del enganche (abanico si comparten vertice)
    por_vertice = {}
    for r, (u, l) in enumerate(patas):
        por_vertice.setdefault(u, []).append(r)
    sig = np_
    etiquetas_t = {}
    for u, rs in por_vertice.items():
        desv = {1: [0], 2: [28, -28], 3: [30, 0, -30]}[len(rs)]
        for r, dv in zip(rs, desv):
            l = patas[r][1]
            a = (float(ang[u]) + dv) * pi / 180
            for k in range(1, l + 1):
                pos[sig] = ((R + paso * k) * cos(a), (R + paso * k) * sin(a))
                sig += 1
            etiquetas_t[ts[r]] = (r, a)

    g = circle((0, 0), R, color="lightgray", linestyle="--", thickness=0.8, zorder=0)
    puentes = set(map(frozenset, G.bridges(labels=False)))
    for u, v in G.edges(labels=False):
        if frozenset((u, v)) in puentes:
            g += line([pos[u], pos[v]], color=NARANJA, thickness=2.6, zorder=1)
        else:
            g += line([pos[u], pos[v]], color="gray", thickness=0.8, zorder=1)
    for v in G.vertices(sort=True):
        if v in ts:
            g += circle(pos[v], 0.11, fill=True, facecolor="white", edgecolor="black", thickness=1.2, zorder=5)
        else:
            g += point(pos[v], size=28, color="black", zorder=4)
    for v, s in etiq.items():
        x, y = pos[v]
        g += text(f"${s}$", (0.66 * x, 0.66 * y), fontsize=11, color="black", zorder=6)
    for t, (r, a) in etiquetas_t.items():
        x, y = pos[t]
        g += text(f"$t_{r + 1}$", (x + 0.42 * cos(a), y + 0.42 * sin(a)), fontsize=12, color="black", zorder=6)
    if titulo is None:
        titulo = f"{caso(n, m)}  $G^*_{{{n},{m}}}$,  $\\ell={terna_justa(b_nm(n, m))}$"
    g += text(titulo, (0, -R - 1.75), fontsize=11, color="black")
    g.set_axes_range(-3.1, 3.1, -3.0, 3.4)
    g.axes(False)
    g.set_aspect_ratio(1)
    return g


# Fila 1: los cuatro casos con b >= 3, como en la figura (ternas (2,1,1) y, en (b'), (2,2,1))
fila_b3 = [dibujar_modelo(n, m) for n, m in [(8, 10), (8, 9), (9, 10), (9, 12)]]
# Fila 2: los mismos casos con b = 2, terna (1,1,0): un terminal queda dentro del nucleo
fila_b2 = [dibujar_modelo(n, m) for n, m in [(7, 12), (7, 11), (7, 10)]]

for n, m in [(8, 10), (8, 9), (9, 10), (9, 12), (7, 12), (7, 11), (7, 10)]:
    G, ts = modelo(n, m)
    b = b_nm(n, m)
    print(f"{caso(n, m):4} G*_{{{n},{m}}}: b={b}, nucleo B_{{{n - b},{m - b}}}, l={terna_justa(b)}, F={F_tupla(G, ts)}")

graphics_array(fila_b3).show(figsize=[16, 4.6])
graphics_array(fila_b2).show(figsize=[12, 4.6])

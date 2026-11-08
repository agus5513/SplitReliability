def DeltaT2H0(b_valor, n_prima, lambda_prima):
    """
        Fórmula exacta para Delta_0 = t2(H_{n,m}^{(0)}) - t2(G_{n,m}),
        válida en el dominio 2 < lambda' <= n' - 2 (mismo dominio que
        NumeroDe2BosquesH0).
    """
    if not (2 < lambda_prima <= n_prima - 2):
        raise ValueError(
            f"DeltaT2H0: fórmula no definida fuera de "
            f"2 < lambda' <= n' - 2 (se obtuvo n'={n_prima}, lambda'={lambda_prima})"
        )
    return (b_valor - 2) * (lambda_prima - 1) * n_prima**(lambda_prima - 2) * (n_prima - 1)**(n_prima - lambda_prima - 1)

def DeltaT2H1(b_valor, n_prima, lambda_prima):
    """
        Fórmula exacta para Delta_1 = t2(H_{n,m}^{(1)}) - t2(G_{n,m}),
        válida en el dominio lambda' <= n' - 4 (mismo dominio que
        NumeroDe2BosquesH1).
    """
    if lambda_prima > n_prima - 4:
        raise ValueError(
            f"DeltaT2H1: fórmula no definida "
            f"para lambda'={lambda_prima} > n' - 4={n_prima - 4}"
        )
    return (
        n_prima**(lambda_prima - 1) * (n_prima - 1)**(n_prima - lambda_prima - 4)
        * (
            (n_prima - 3) * (lambda_prima * (n_prima - 1) * b_valor + n_prima + lambda_prima)
            - 2 * lambda_prima * (n_prima - 1) * (n_prima - 2)
        )
    )

def DeltaT2H2(b_valor, n_prima, lambda_prima):
    """
        Fórmula exacta para Delta_2 = t2(H_{n,m}^{(2)}) - t2(G_{n,m}),
        válida en el dominio 2 <= lambda' <= n' - 2 (mismo dominio que
        NumeroDe2BosquesH2).
    """
    if not (2 <= lambda_prima <= n_prima - 2):
        raise ValueError(
            f"DeltaT2H2: fórmula no definida fuera de "
            f"2 <= lambda' <= n' - 2 (se obtuvo n'={n_prima}, lambda'={lambda_prima})"
        )
    return (
        n_prima**(lambda_prima - 2) * (n_prima - 1)**(n_prima - lambda_prima - 3)
        * (
            n_prima * (n_prima - lambda_prima - 2)
            + lambda_prima * (n_prima - 2) * (n_prima - 1) * (b_valor - 2)
        )
    )

def DeltaT2(k, b_valor, n_prima, lambda_prima):
    """
        Función wrapper que evalúa Delta_k (k = 0, 1, 2) delegando
        en DeltaT2H0, DeltaT2H1 o DeltaT2H2 según corresponda.
    """
    if k == 0:
        return DeltaT2H0(b_valor, n_prima, lambda_prima)
    elif k == 1:
        return DeltaT2H1(b_valor, n_prima, lambda_prima)
    elif k == 2:
        return DeltaT2H2(b_valor, n_prima, lambda_prima)
    raise ValueError(f"k={k} no es un índice válido (debe ser 0, 1 o 2)")

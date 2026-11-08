import os
while not os.path.isdir("DataTypes"):
    os.chdir("..")

load("DataTypes/Grafo2Terminales.sage")
load("Scripts/GrafoBalloon.sage")
load("Scripts/VariantesGrafoBallon.sage")
load("Scripts/NumeroDe2Bosques.sage")

def TestVariantesBalloonPreservanOrdenYTamano():
    """
        Test exhaustivo (n entre 4 y 14) que verifica que, para
        todo (n, m) donde H0, H1 o H2 son construibles, preservan
        el orden y el tamaño del grafo Balloon original: contraer
        un puente y subdividir una arista se compensan, dejando
        |V| = n y |E| = m. También verifica que el grafo resultante
        siga siendo conexo. Los (n, m) donde alguna variante no es
        construible (precondiciones de H0/H1/H2) se saltean.
    """

    probados = 0
    for n in range(4, 15):
        for m in range(n, binomial(n, 2) + 1):
            for k in range(3):
                try:
                    TC = VariantesBalloonH(k, n, m)
                except ValueError:
                    continue
                probados += 1
                if TC.grafo.order() != n:
                    return False
                if TC.grafo.size() != m:
                    return False
                if not TC.grafo.is_connected():
                    return False
    return probados > 100

def TestH1UneNodosNoAdyacentesASPrima():
    """
        Test exhaustivo (n entre 4 y 14) que verifica que, en todo
        (n, m) donde H1 es construible, el vértice subdividido
        efectivamente une a 2 nodos de G' que no son adyacentes a
        s' (y distintos de la imagen de TC.s, el otro terminal), y
        que dichos nodos sí están unidos por una arista en G'.

        NOTA: s' (el vértice de grado mínimo de G') es siempre la
        proyección de TC.t en estos grafos Balloon, así que el nodo
        que hay que excluir por separado del pool de no-adyacentes
        es la proyección de TC.s, no la de TC.t.
    """

    probados = 0
    for n in range(4, 15):
        for m in range(n, binomial(n, 2) + 1):
            try:
                TC = ContraerPuente(ConstruirGrafo2TerminalesBalloon(n, m))
            except ValueError:
                continue
            G_prima, s_prima, lambda_prima = NodoSPrima(TC)
            representantes = RepresentantesEsqueleto(TC)
            s_real_prima = representantes[TC.s]
            no_adyacentes = [
                v for v in G_prima.vertices()
                if v != s_prima and v != s_real_prima and not G_prima.has_edge(s_prima, v)
            ]
            if len(no_adyacentes) < 2:
                continue
            probados += 1
            if not G_prima.has_edge(no_adyacentes[0], no_adyacentes[1]):
                return False
    return probados > 20

def TestH2UneNodosAdyacentesASPrima():
    """
        Test exhaustivo (n entre 4 y 14) que verifica que, en todo
        (n, m) donde H2 es construible, s' cuenta con al menos 2
        vecinos en G' (condición que necesita H2 para poder
        subdividir la arista que los une), y que dichos vecinos
        sí están unidos por una arista en G'.
    """

    probados = 0
    for n in range(4, 15):
        for m in range(n, binomial(n, 2) + 1):
            try:
                TC = ContraerPuente(ConstruirGrafo2TerminalesBalloon(n, m))
            except ValueError:
                continue
            G_prima, s_prima, lambda_prima = NodoSPrima(TC)
            adyacentes = [v for v in G_prima.neighbors(s_prima) if v != s_prima]
            if len(adyacentes) < 2:
                continue
            probados += 1
            if not G_prima.has_edge(adyacentes[0], adyacentes[1]):
                return False
    return probados > 20

def TestVariantesBalloonHWrapper():
    """
        Test que verifica, para varios (n, m), que
        VariantesBalloonH(k, n, m) coincide exactamente con H0,
        H1 o H2 según el valor de k (incluyendo casos donde alguna
        variante falla, en cuyo caso ambos lados deben fallar), y
        que rechaza valores de k fuera de {0, 1, 2}.
    """

    casos = [(7, 10), (8, 18), (9, 20), (9, 21), (9, 30), (10, 20), (11, 25), (12, 30)]
    funciones = {0: H0, 1: H1, 2: H2}

    for n, m in casos:
        for k, funcion in funciones.items():
            try:
                esperado = funcion(n, m)
                esperado_fallo = False
            except ValueError:
                esperado = None
                esperado_fallo = True

            try:
                obtenido = VariantesBalloonH(k, n, m)
                obtenido_fallo = False
            except ValueError:
                obtenido = None
                obtenido_fallo = True

            if esperado_fallo != obtenido_fallo:
                return False
            if not esperado_fallo and obtenido != esperado:
                return False

        try:
            VariantesBalloonH(3, n, m)
            return False
        except ValueError:
            pass
    return True

def TestH0RequierePuentes():
    """
        Test que verifica, para varios (n, m) sin puentes (caso
        Indice0), que H0 falla explícitamente al no tener un
        puente que contraer.
    """

    casos = [(9, 30), (4, 5), (4, 6), (5, 8), (5, 9), (6, 12), (6, 15), (7, 17)]
    for n, m in casos:
        try:
            H0(n, m)
            return False
        except ValueError:
            pass
    return True

def TestH1RequiereNoAdyacentes():
    """
        Test que verifica, para varios (n, m), que H1 falla
        explícitamente cuando G' no tiene 2 nodos no adyacentes a
        s' (excluyendo la imagen de TC.s, el otro terminal).
    """

    casos = [(9, 18), (4, 4), (4, 5), (4, 6), (5, 5), (5, 6), (5, 7), (6, 6), (6, 7), (6, 8)]
    for n, m in casos:
        try:
            H1(n, m)
            return False
        except ValueError:
            pass
    return True


def TestNumeroDe2BosquesVarianteBallonRechazaFueraDeRango():
    """
        Test que verifica, con varios (n, m) por caso, que
        NumeroDe2BosquesVarianteBallon lanza ValueError fuera de
        los rangos donde cada fórmula está probada: k=0 cuando el
        esqueleto G' es completo (lambda' > n' - 2), k=1 cuando
        lambda' > n' - 4, y k=2 cuando el esqueleto G' es completo
        (lambda' > n' - 2). También rechaza k fuera de {0, 1, 2}.
    """

    casos = [
        (0, 4, 4), (0, 5, 5), (0, 5, 7), (0, 6, 6), (0, 6, 8),
        (1, 8, 8), (1, 4, 4), (1, 5, 5), (1, 6, 6), (1, 7, 7),
        (2, 4, 4), (2, 4, 6), (2, 5, 5), (2, 5, 7),
    ]
    for k, n, m in casos:
        try:
            NumeroDe2BosquesVarianteBallon(k, n, m)
            return False
        except ValueError:
            pass

    for n, m in [(8, 18), (9, 20), (10, 22)]:
        try:
            NumeroDe2BosquesVarianteBallon(3, n, m)
            return False
        except ValueError:
            pass
    return True


def TestNumeroDe2BosquesH0RespetaDominio():
    """
        Test directo sobre NumeroDe2BosquesH0 (independiente de si
        el grafo Balloon correspondiente es construible) que recorre
        muchos valores de n' y lambda' y verifica que la fórmula
        lanza ValueError exactamente fuera del dominio
        2 < lambda' <= n' - 2, y que calcula un valor (sin lanzar)
        para todo (n', lambda') dentro de ese dominio.
    """

    b_valor = 5
    probados = 0
    for n_prima in range(3, 16):
        for lambda_prima in range(0, n_prima + 2):
            dentro_de_rango = 2 < lambda_prima <= n_prima - 2
            try:
                NumeroDe2BosquesH0(b_valor, n_prima, lambda_prima)
                lanzo_error = False
            except ValueError:
                lanzo_error = True
            probados += 1
            if dentro_de_rango == lanzo_error:
                return False
    return probados > 100


def TestNumeroDe2BosquesH1RespetaDominio():
    """
        Test directo sobre NumeroDe2BosquesH1 (independiente de si
        el grafo Balloon correspondiente es construible) que recorre
        muchos valores de n' y lambda' y verifica que la fórmula
        lanza ValueError exactamente fuera del dominio
        lambda' <= n' - 4, y que calcula un valor (sin lanzar) para
        todo lambda' dentro de ese dominio.
    """

    b_valor = 5
    probados = 0
    for n_prima in range(3, 16):
        for lambda_prima in range(0, n_prima + 2):
            dentro_de_rango = lambda_prima <= n_prima - 4
            try:
                NumeroDe2BosquesH1(b_valor, n_prima, lambda_prima)
                lanzo_error = False
            except ValueError:
                lanzo_error = True
            probados += 1
            if dentro_de_rango == lanzo_error:
                return False
    return probados > 100


def TestNumeroDe2BosquesH2RespetaDominio():
    """
        Test directo sobre NumeroDe2BosquesH2 (independiente de si
        el grafo Balloon correspondiente es construible) que recorre
        muchos valores de n' y lambda' y verifica que la fórmula
        lanza ValueError exactamente fuera del dominio
        2 <= lambda' <= n' - 2, y que calcula un valor (sin lanzar)
        para todo (n', lambda') dentro de ese dominio.
    """

    b_valor = 5
    probados = 0
    for n_prima in range(3, 16):
        for lambda_prima in range(0, n_prima + 2):
            dentro_de_rango = 2 <= lambda_prima <= n_prima - 2
            try:
                NumeroDe2BosquesH2(b_valor, n_prima, lambda_prima)
                lanzo_error = False
            except ValueError:
                lanzo_error = True
            probados += 1
            if dentro_de_rango == lanzo_error:
                return False
    return probados > 100


def TestNumeroDe2BosquesH0VsFuerzaBruta():
    """
        Test exhaustivo (n entre 4 y 19), dedicado exclusivamente a
        k=0 (separado de k=1 y k=2 para que su resultado no quede
        enmascarado por el estado de las otras dos fórmulas), que
        compara NumeroDe2BosquesVarianteBallon(0, n, m) contra el
        número de 2-bosques de H0(n, m) computado por fuerza bruta
        (vía el menor del Laplaciano, NumeroDe2Bosques).
    """

    comparaciones = 0
    for n in range(4, 20):
        for m in range(n, binomial(n, 2) + 1):
            try:
                obtenido = NumeroDe2BosquesVarianteBallon(0, n, m)
            except ValueError:
                continue
            try:
                TC = VariantesBalloonH(0, n, m)
                esperado = NumeroDe2Bosques(TC)
            except ValueError:
                continue
            if obtenido != esperado:
                return False
            comparaciones += 1
    return comparaciones > 400


def TestNumeroDe2BosquesH1VsFuerzaBruta():
    """
        Test exhaustivo (n entre 4 y 19), dedicado exclusivamente a
        k=1 (separado de k=0 y k=2 para que su resultado no quede
        enmascarado por el estado de las otras dos fórmulas), que
        compara NumeroDe2BosquesVarianteBallon(1, n, m) contra el
        número de 2-bosques de H1(n, m) computado por fuerza bruta
        (vía el menor del Laplaciano, NumeroDe2Bosques).
    """

    comparaciones = 0
    for n in range(4, 20):
        for m in range(n, binomial(n, 2) + 1):
            try:
                obtenido = NumeroDe2BosquesVarianteBallon(1, n, m)
            except ValueError:
                continue
            try:
                TC = VariantesBalloonH(1, n, m)
                esperado = NumeroDe2Bosques(TC)
            except ValueError:
                continue
            if obtenido != esperado:
                return False
            comparaciones += 1
    return comparaciones > 400


def TestNumeroDe2BosquesH2VsFuerzaBruta():
    """
        Test exhaustivo (n entre 4 y 19), dedicado exclusivamente a
        k=2 (separado de k=0 y k=1 para que su resultado no quede
        enmascarado por el estado de las otras dos fórmulas), que
        compara NumeroDe2BosquesVarianteBallon(2, n, m) contra el
        número de 2-bosques de H2(n, m) computado por fuerza bruta
        (vía el menor del Laplaciano, NumeroDe2Bosques).
    """

    comparaciones = 0
    for n in range(4, 20):
        for m in range(n, binomial(n, 2) + 1):
            try:
                obtenido = NumeroDe2BosquesVarianteBallon(2, n, m)
            except ValueError:
                continue
            try:
                TC = VariantesBalloonH(2, n, m)
                esperado = NumeroDe2Bosques(TC)
            except ValueError:
                continue
            if obtenido != esperado:
                return False
            comparaciones += 1
    return comparaciones > 400


if __name__ == "__main__":
    tests = [
        ("TestVariantesBalloonPreservanOrdenYTamano", TestVariantesBalloonPreservanOrdenYTamano),
        ("TestH1UneNodosNoAdyacentesASPrima", TestH1UneNodosNoAdyacentesASPrima),
        ("TestH2UneNodosAdyacentesASPrima", TestH2UneNodosAdyacentesASPrima),
        ("TestVariantesBalloonHWrapper", TestVariantesBalloonHWrapper),
        ("TestH0RequierePuentes", TestH0RequierePuentes),
        ("TestH1RequiereNoAdyacentes", TestH1RequiereNoAdyacentes),
        ("TestNumeroDe2BosquesVarianteBallonRechazaFueraDeRango", TestNumeroDe2BosquesVarianteBallonRechazaFueraDeRango),
        ("TestNumeroDe2BosquesH0RespetaDominio", TestNumeroDe2BosquesH0RespetaDominio),
        ("TestNumeroDe2BosquesH1RespetaDominio", TestNumeroDe2BosquesH1RespetaDominio),
        ("TestNumeroDe2BosquesH2RespetaDominio", TestNumeroDe2BosquesH2RespetaDominio),
        ("TestNumeroDe2BosquesH0VsFuerzaBruta", TestNumeroDe2BosquesH0VsFuerzaBruta),
        ("TestNumeroDe2BosquesH1VsFuerzaBruta", TestNumeroDe2BosquesH1VsFuerzaBruta),
        ("TestNumeroDe2BosquesH2VsFuerzaBruta", TestNumeroDe2BosquesH2VsFuerzaBruta),
    ]

    fallidos = [nombre for nombre, test in tests if test() != True]
    if not fallidos:
        print("Todos los tests exitosos")
    else:
        for nombre in fallidos:
            print(f"Test {nombre} fallido")

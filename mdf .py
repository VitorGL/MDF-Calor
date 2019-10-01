from threading import *
from queue import Queue
from time import *

q = None
fourier = 0
w_lock = Lock()
solido = []
solido2 = []


def main():
    h = 1
    ap = 0
    tam = 10
    ans = None
    coef_cond = 1.0

    while True:

        ans = int(input("Precisao:\n"
            "(1) - baixa\n"
            "(2) - media\n"
            "(3) - alta\n"))

        if ans == 1:
            ap = 2.0
            break

        elif ans == 2:

            ap = 4.0
            break

        elif ans == 3:

            ap = 6.0
            break

        else:
            print("resposta inesperada")

    tempo = h**2/ap*coef_cond+1

    while tempo > h**2/ap*coef_cond:

        tempo = float(input("Digite um valor para tempo, tal que 'tempo < {}'".format(h**2/ap*coef_cond)))

    print("Considere o valor da condutividade termica em cm^2/s")

    calor_mdf_thread(h, tempo, tam, coef_cond)


def cria_matriz(x, y, z, value):

    matriz = []

    for i in range(x):
        matriz.append([])
        for j in range(y):
            matriz[i].append([])
            for k in range(z):
                matriz[i][j].append(value)

    return matriz


def mod_temp_plano(m, pos, temp):

    t = len(m[pos])/2

    if t % 2 == 0:
        nj = t/2
    else:
        nj = (t-1)/2

    t = len(m[pos][0])/2

    if t % 2 == 0:
        nk = t / 2
    else:
        nk = (t - 1) / 2

    for j in range(len(m[pos])):
        for k in range(len(m[pos][j])):
            if nj-1 < j < len(m[pos][j])-nj and nk-1 < k < len(m[pos][j])-nk:
                m[pos][j][k] = temp

    return m


def calor_mdf(h, tempo, dimensao, alfa):
    mod_temp = float(input("Digite a temperatura que sera aplicada a uma area do cubo:"))

    inicio = time()
    global solido
    global solido2
    global fourier

    d2 = dimensao+2
    temp = 35
    mod_temp = 0

    # print("Digite a temperatura base do cubo:\n")
    # scanf("%lf", &temp)

    solido = mod_temp_plano(cria_matriz(d2, d2, d2, 35), 0, mod_temp)

    # solido2 = solido[:]
    solido2 = [[solido[i][j][:] for j in range(len(solido[i]))] for i in range(len(solido))]

    on = True
    c_vizinhas_som = 0
    Z = 0
    Z2 = 0
    valor = 0

    fourier = alfa**2 * (tempo/h**2)

    while on:
        # solido2 = solido[:]
        solido2 = [[solido[i][j][:] for j in range(len(solido[i]))] for i in range(len(solido))]

        for i in range(1, d2-1):
            for j in range(1, d2-1):
                for k in range(1, d2-1):

                    if 0 < i <= dimensao and 0 < j <= dimensao and 0 < k <= dimensao:

                        c_vizinhas_som = (solido[i][j+1][k]
                                       + solido[i][j-1][k]
                                       + solido[i-1][j][k]
                                       + solido[i+1][j][k]
                                       + solido[i][j][k+1]
                                       + solido[i][j][k-1]
                        )

                        solido2[i][j][k] = solido[i][j][k] + fourier * (c_vizinhas_som - (6 * solido[i][j][k]))  # Equação

        # print_matriz(solido2)

        valor = 0

        # for (int i = 0 i < d2 i++)
        #     for (int j = 0 j < d2 j++)
        #         for (int k = 0 k < d2 k++)
        #
        #             valor += abs(solido[i][j][k] - solido2[i][j][k])
        #             # print("valor = %lf = |%lf - %lf|\n\n", valor, solido[i][j][k], solido2[i][j][k])
        #             # sleep(1)
        #

        for i in range(dimensao):
            for j in range(dimensao):
                for k in range(dimensao):

                    valor += abs(solido[i][j][k] - solido2[i][j][k])
                    # print("valor = %lf = |%lf - %lf|\n\n", valor, solido[i][j][k], solido2[i][j][k])
                    # sleep(1)

        Z = valor / (d2 * d2 * d2)

        print("Z = {}\n".format(Z))
        # sleep(1)

        if Z <= 0:  # || Z == Z2

            print("A condicaoo de estabilidade foi atingida com valor ({} < 0)".format(Z))
            on = False

        else:
            Z2 = Z

        # solido = solido2[:]
        solido = [[solido2[i][j][:] for j in range(len(solido2[i]))] for i in range(len(solido2))]
    print("Tempo=", time() - inicio)


def print_matriz(m):

    for i in range(len(m)):
        for j in range(len(m[i])):
            for k in range(len(m[i][j])):
                # print("M(%d)(%d)(%d) = %lf  ", i, j, k, m[i][j][k])
                print("{0:.8f}  ".format(m[i][j][k]), end="")
            print()

        print("\n")


def calor_mdf_thread(h, tempo, dimensao, alfa):
    mod_temp = float(input("Digite a temperatura que sera aplicada a uma area do cubo:"))

    inicio = time()
    global solido
    global solido2
    global fourier

    d2 = dimensao + 2
    temp = 35

    # print("Digite a temperatura base do cubo:\n")
    # scanf("%lf", &temp)

    solido = mod_temp_plano(cria_matriz(d2, d2, d2, 35), 1, mod_temp)

    solido2 = [[solido[i][j][:] for j in range(len(solido[i]))] for i in range(len(solido))]

    on = True
    c_vizinhas_som = 0
    Z = 0
    Z2 = 0
    valor = 0

    fourier = alfa ** 2 * (tempo / h ** 2)

    ini_thread = time()

    global q
    q = Queue()

    for x in range(3):
        t = Thread(target=threader)
        t.daemon = True
        t.start()

    qntx = 2
    qnty = 2
    qntz = 2

    subx = int((len(solido)-2) / qntx)
    suby = int((len(solido[0])-2) / qnty)
    subz = int((len(solido[0][0])-2) / qntz)

    intervalos = []

    for i in range(qntx):
        cx = (i * subx)+1
        fx = ((i * subx) + subx)
        # print(cx, fx)
        for j in range(qnty):
            cy = (j * suby)+1
            fy = ((j * suby) + suby)
            # print(cy, fy)
            for k in range(qntz):
                cz = (k * subz)+1
                fz = ((k * subz) + subz)
                # print(cz, fz)

                # matriz = [[solido[i][j][cz:fz] for j in range(cy, fy)] for i in range(cx, fx)]
                # print(matriz)
                intervalos.append([[cx, fx], [cy, fy], [cz, fz]])

    tempo_thread = time() - ini_thread

    while on:
        solido2 = [[solido[i][j][:] for j in range(len(solido[i]))] for i in range(len(solido))]

        for intervalo in intervalos:
            q.put(intervalo)

        q.join()
        # print_matriz(solido2)

        valor = 0

        for i in range(1, dimensao):
            for j in range(1, dimensao):
                for k in range(1, dimensao):
                    valor += abs(solido[i][j][k] - solido2[i][j][k])

        Z = valor / (d2 * d2 * d2)

        print("Z = {}\n".format(Z))

        if Z <= 0:  # || Z == Z2

            print("A condicaoo de estabilidade foi atingida com valor ({} < 0)".format(Z))
            on = False

        else:
            Z2 = Z

        solido = [[solido2[i][j][:] for j in range(len(solido2[i]))] for i in range(len(solido2))]

    print("Tempo=", time() - inicio)
    print("Tempo de inicialização a mais das threads=", tempo_thread)


def calculo_do_ponto(coords):
    global solido
    global solido2
    global fourier
    cx, fx = coords[0]
    cy, fy = coords[1]
    cz, fz = coords[2]

    '''Calculo da função'''
    for i in range(len(solido)):
        for j in range(len(solido[i])):
            for k in range(len(solido[i][j])):
                if cx <= i <= fx and cy <= j <= fy and cz <= k <= fz:
                    c_vizinhas_som = (solido[i][j + 1][k]
                                      + solido[i][j - 1][k]
                                      + solido[i - 1][j][k]
                                      + solido[i + 1][j][k]
                                      + solido[i][j][k + 1]
                                      + solido[i][j][k - 1]
                                      )

                    # colocar lock aqui
                    w_lock.acquire()
                    solido2[i][j][k] = solido[i][j][k] + fourier * (c_vizinhas_som - (6 * solido[i][j][k]))
                    w_lock.release()


def threader():
    global q
    while True:
        coords = q.get()
        calculo_do_ponto(coords)
        q.task_done()


main()
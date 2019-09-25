
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

    calorMDF(h, tempo, tam, coef_cond)


def criaMatriz(x, y, z, value):

    matriz = []

    for i in range(x):
        matriz.append([])
        for j in range(y):
            matriz[i].append([])
            for k in range(z):
                matriz[i][j].append(value)
            
    return matriz


def modTempPlane(m, x, y, z, pos, temp):

    t = y/2

    if t % 2 == 0:
        nj = t/2
    else:
        nj = (t-1)/2

    t = z / 2

    if t % 2 == 0:
        nk = t / 2
    else:
        nk = (t - 1) / 2

    for j in range(y):
        for k in range(z):
            if nj-1 < j < z-nj and nk-1 < k < z-nk:
                m[pos][j][k] = temp
        
    return m


def calorMDF(h, tempo, dimensao, alfa):

    d2 = dimensao+2
    temp = 35
    mod_temp = 0

    # print("Digite a temperatura base do cubo:\n")
    # scanf("%lf", &temp)

    mod_temp = float(input("Digite a temperatura que sera aplicada a uma area do cubo:"))

    u = modTempPlane(criaMatriz(d2, d2, d2, 35), d2, d2, d2, 0, mod_temp)

    u2 = u[:]

    on = True
    c_vizinhas_som = 0
    Z = 0
    Z2 = 0
    valor = 0

    fourier = alfa**2 * (tempo/h**2)

    while on:
        u2 = u[:]

        for i in range(dimensao):
            for j in range(dimensao):
                for k in range(dimensao):
                
                    if 0 < i <= dimensao and 0 < j <= dimensao and 0 < k <= dimensao:
                    
                        c_vizinhas_som = (u[i][j+1][k]
                                       + u[i][j-1][k]
                                       + u[i-1][j][k]
                                       + u[i+1][j][k]
                                       + u[i][j][k+1]
                                       + u[i][j][k-1]
                        )

                        u2[i][j][k] = u[i][j][k] + fourier * (c_vizinhas_som - (6 * u[i][j][k]))  # Equação

        print_matriz(u2)

        valor = 0

        # for (int i = 0 i < d2 i++)
        #     for (int j = 0 j < d2 j++)
        #         for (int k = 0 k < d2 k++)
        #         
        #             valor += abs(u[i][j][k] - u2[i][j][k])
        #             # print("valor = %lf = |%lf - %lf|\n\n", valor, u[i][j][k], u2[i][j][k])
        #             # sleep(1)
        #         

        for i in range(dimensao):
            for j in range(dimensao):
                for k in range(dimensao):
                
                    valor += abs(u[i][j][k] - u2[i][j][k])
                    # print("valor = %lf = |%lf - %lf|\n\n", valor, u[i][j][k], u2[i][j][k])
                    # sleep(1)

        Z = valor / (d2 * d2 * d2)

        print("Z = {}\n".format(Z))
        # sleep(1)

        if Z <= 0:  # || Z == Z2
        
            print("A condicaoo de estabilidade foi atingida com valor ({} < 0)".format(Z))
            on = False
        
        else:
            Z2 = Z

        u = u2[:]

    # FILE *arq

    # if ((arq = fopen("log.txt", "w")) == NULL)
    #     print("Erro ao abrir arquivo de log\n")

    # for (int i = 1 i <= dimensao i++)
    
    #     for (int j = 1 j <= dimensao j++)
        
    #         for (int k = 1 k <= dimensao k++)
    #             fprint(arq, "%lf  ", u[i][j][k])
    #         fprint(arq, "\n")
        
    #     fprint(arq, "\n\n")

    # fclose(arq)


def print_matriz(m):

    for i in range(len(m)):
        for j in range(len(m[i])):
            for k in range(len(m[i][j])):
                # print("M(%d)(%d)(%d) = %lf  ", i, j, k, m[i][j][k])
                print("{0:f}  ".format(m[i][j][k]), end="")
            print()
        
        print("\n")


def calorMDFThread(h, tempo, dimensao, alfa):

    d2 = dimensao + 2
    temp = 35

    # print("Digite a temperatura base do cubo:\n")
    # scanf("%lf", &temp)

    mod_temp = float(input("Digite a temperatura que sera aplicada a uma area do cubo:"))

    u = modTempPlane(criaMatriz(d2, d2, d2, 35), d2, d2, d2, 0, mod_temp)

    u2 = u[:]

    on = True
    c_vizinhas_som = 0
    Z = 0
    Z2 = 0
    valor = 0

    fourier = alfa ** 2 * (tempo / h ** 2)

    global q
    q = Queue()

    for x in range(min(dimensao, 100)):
        t = Thread(target=threader)
        t.daemon = True
        t.start()

    while on:
        u2 = u[:]

        for i in range(dimensao):
            q.put(u[i])

        q.join()

        print_matriz(u2)

        valor = 0

        for i in range(1, dimensao):
            for j in range(1, dimensao):
                for k in range(1, dimensao):
                    valor += abs(u[i][j][k] - u2[i][j][k])

        Z = valor / (d2 * d2 * d2)

        print("Z = {}\n".format(Z))

        if Z <= 0:  # || Z == Z2

            print("A condicaoo de estabilidade foi atingida com valor ({} < 0)".format(Z))
            on = False

        else:
            Z2 = Z

        u = u2[:]

def calculo_do_ponto(i):
    '''Calculo da função'''
    for j in range(dimensao):
        for k in range(dimensao):
            if 0 < i <= dimensao and 0 < j <= dimensao and 0 < k <= dimensao:
                c_vizinhas_som = (u[i][j + 1][k]
                                  + u[i][j - 1][k]
                                  + u[i - 1][j][k]
                                  + u[i + 1][j][k]
                                  + u[i][j][k + 1]
                                  + u[i][j][k - 1]
                                  )

                u2[i][j][k] = u[i][j][k] + fourier * (c_vizinhas_som - (6 * u[i][j][k]))




def threader():
    global q
    while True:
        x = q.get()
        calculo_do_ponto(x)
        q.task_done()

# exemplo
def runPS():
    global q
    q = Queue()

    for x in range(100):
        t = Thread(target=threader)
        t.daemon = True
        t.start()

    tam = len(m) + len(m[]) + len(m[][])

    for worker in range(4000, 5000):
        q.put(worker)

    q.join()

main()
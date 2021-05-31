'''
Este script define as funções que criam os ângulos dos N polinômios
de cada junta.

Este arquivo pode ser importado como um módulo utilizando:
from src.utils.numerical import create_angles
'''

import random
import numpy as np

def create_angles(theta_i, theta_f, final_time, steps):
    '''
    Função para definição dos múltiplos polinômios presentes
    para cada junta.

    Parâmetros
    ----------
    theta_i: list
        Lista contendo posições angulares iniciais de todas juntas.
    theta_f: list
        Lista contendo posições angulares finais de todas juntas.
    final_time : int
        Tempo de duração da simulação.
    steps : int
        Inteiro que representa o número de polinômios que serão
        criados.

    Parâmetros
    ----------
    number : int
        Inteiro que representa o número de pontos do polinômio.
    '''
    theta, omega, time = [theta_i], [0], [0]
    for i in range(steps-2):
        theta.append((theta_f - theta_i)*random.random() + theta_i)
        omega.append(5*random.random())
        time.append(final_time*random.random())
    omega.append(0)
    theta.append(theta_f)
    time.append(final_time)
    time.sort()
    if theta_i < theta_f:
        theta.sort()
    else:
        theta.sort(reverse=True)
    return time, theta, omega


        
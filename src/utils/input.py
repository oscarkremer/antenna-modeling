'''
Este script define as funções de input, tanto de variáveis
no espaço de juntas quanto espaço cartesiano.

Este arquivo pode ser importado como um módulo utilizando:
from src.utils.input import input_angles, input_cartesian
'''
import numpy as np

def input_angles():
    '''
    Função para input de ângulos do robô para cálculo da cinemática
    inversa. Os angulos das juntas são inseridos em graus, convertidos
    para radianos e, por fim, retornados. 

    Retorna
    -------
    angles
        Lista com os ângulos das juntas do robô.
    '''
    angles = np.zeros(5)
    for i in range(5):
        angles[i] = input('-- Theta {}: '.format(i+1))
        angles[i] = np.pi*angles[i]/180
    return angles

def input_cartesian():
    '''
    Função para input das variaveis de posição e orientação do atuador 
    final do robo no plano cartesiano. Os ângulos alpha, beta e gamma
    são tratados em radianos mas entrados em graus. 

    Retorna
    -------
    pos
        Booleano que será verdadeiro caso o par posição/orientação
        inseridas não seja alcançável no espaço de tarefa.
    angles
        Lista com os ângulos da orientação do robô (alpha, beta, gamma).
    '''
    pos = np.zeros(3)
    angles = np.zeros(3)
    pos[0] = input('-- x (cm): ')
    pos[1] = input('-- y (cm): ')
    pos[2] = input('-- z (cm): ')
    angles[0] = np.pi*(float(input('-- alpha: ')))/180
    angles[1] = np.pi*(float(input('-- beta: ')))/180
    angles[2] = np.pi*(float(input('-- gamma: ')))/180
    return pos, angles
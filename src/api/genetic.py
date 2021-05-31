'''
Este script permite o treinamento do algoritmo genético para 
otimização do robô mentor, assim como publicado no trabalho:  
A Genetic Approach for Trajectory Optimization Applied to a 
Didactic Robot. O. S. Kremer; M. A. B. Cunha; F. S. Moraes; 
S. S. Schiavon. 2019 Latin American Robotics Symposium (LARS), 
2019 Brazilian Symposium on Robotics (SBR).

O script pede a inserção da posição inicial e final do atuador 
final e ângulos finais e iniciais para orientação do mesmo.

Com o ambiente mentor ativado no conda este script pode ser 
executado com o comando make genetic.
'''
import argparse
import numpy as np
from src.models import Population
from src.mentor import Mentor
from src.utils.input import input_cartesian


def calculate_thetas(pos, angles):
    '''
    Função para calcular ângulos das juntas a partir
    de ângulos de orientação e posição.

    Parâmetros
    ----------
    pos : list
        Lista com as posições (x, y, z) no espaço cartesiano.
    angles: list
        Lista com os ângulos que representam a orientação XYZ 
        do atuador final.
    Retorna
    -------
    error
        Booleano que será verdadeiro caso a posição/orientação
        inseridas não sejam alcançáveis no espaço de tarefa.
    theta
        Lista com os ângulos do robô para a posição.
    '''
    robot = Mentor()
    rot  = robot.get_orientation(angles[0], angles[1], angles[2])
    matrix_G0 = [[rot[0][0], rot[0][1], rot[0][2], pos[0]],
    [rot[1][0], rot[1][1], rot[1][2], pos[1]],
    [rot[2][0], rot[2][1], rot[2][2], pos[2]],
    [0, 0, 0, 1]]
    matrix_5G = [[1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, -5],
    [0, 0, 0, 1]]
    matrix = np.matmul(matrix_G0, matrix_5G)
    positions = [matrix[0][3], matrix[1][3], matrix[2][3]]
    error, theta = robot.get_angles(positions, rot)
    return error, theta


def enter_position():
    '''
    Função para entrar com os ângulo e posição e testar 
    se são possíveis, mantendo a inserção em um loop até
    que condições possível sejam inseridas

    Parâmetros
    ----------
    Não há parâmetros de entrada nesta função

    Retorna
    -------
    theta
        Lista com os ângulos das juntas encontrados na 
        cinemática inversa.
    pos
        Lista contendo as coordenadas cartesianas que definem
        a posição inserida.
    '''
    error = True
    while error:
        pos, angles = input_cartesian()    
        error, theta = calculate_thetas(pos, angles)
        if error:
            print('Erro !!! \n Posição e/ou orientação não é atingível pelo atuador !!! \n Por favor teste outros valores!')
    return theta, pos


def main(steps, time, generations, mode, population, cross_over, mutation):
    '''
    Função principal que inclui inserção de posição e orientação desejada, 
    instanciação da população para otimização e execução do algoritmo de otimização.

    Parâmetros
    ----------
    steps: int
        Número de sub-polinômios que definem as curvas de trajetória das juntas.
    time: float
        Tempo que irá durar o movimento executado pelo robô Mentor.
    generations: int
        Número de gerações que serão utilizados.
    mode: str
        Modo de adaptação das probabilidades.
    population: int
        Tamanho da população para definir etapa de seleção.
    cross_over: float
        Probabilidade de ocorrer cross-over entre dois elementos distintos.
    mutation: float
        Probabilidade de ocorrer mutação nos elementos.
    '''
    theta_i, pos_i = enter_position()
    theta_f, pos_f = enter_position()
    optimized = Population(population, population, generations, cross_over, mutation, mode=mode)
    population = optimized.initialization(theta_i, theta_f, time, steps)
    optimized.generation(population, theta_i, theta_f, time, steps)


if __name__ == "__main__":
    parser = argparse.ArgumentParser('Inicializando Algoritmo de Otimização')
    parser.add_argument('--mode', default='normal', type=str, 
        choices=['normal', 'adaptive'])
    parser.add_argument('--steps', default=3, type=int)
    parser.add_argument('--time', default=10, type=float)
    parser.add_argument('--generations', default=15, type=int)
    parser.add_argument('--population', default=15, type=int)
    parser.add_argument('--cross-over', default=0.9, type=float)
    parser.add_argument('--mutation', default=0.5, type=float)
    args = parser.parse_args()
    main(args.steps, args.time, args.generations, args.mode, args.population, args.cross_over, args.mutation)

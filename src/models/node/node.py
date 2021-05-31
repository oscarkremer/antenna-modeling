'''
Este script define a classe de cada elemento da população 

Este arquivo pode ser importado como um módulo utilizando:
from src.models.node import Node

Esta importação é feita apenas dentro da pasta models. Para
alterar esta condição modifique o __init__ que encontra-se 
junto a este arquivo.
'''

import numpy as np
from src.mentor import Mentor
from src.polinomy import Polinomy
from src.utils.numerical import create_angles
from src.utils.constants import MAXIMUM_VELOCITY, POINTS

class Node:
    '''
    Classe de cada elemento da população. A classe será definida
    por elementos que obedeçam as leis de cinemática do robô mentor,
    além de limitações (constraints) de velocidade e posição.

    Atributos
    ----------
    joint : list
        List contendo tempo, angulos e velocidades discretizadas.
    dist: int
        Métrica de distância percorrida do elemento.
    points: int
        Coordenadas cartesianos do elemento ao longo do tempo.
    angle: Numpy Array
        Vetor bidimensional dos ângulos de todas as juntas ao
        longo do tempo.
    constraint: bool
        Tag para identificar violação de constraints por parte do
        elemento. Velocidade máxima permitida definida em src/utils/constants.py

    Métodos
    -------
    find_points(self, parents, theta_i, theta_f, time, steps):
        Efetua a operação de cross-over com base em dois pais de um membro.

    test_velocity(self, member, theta_i, theta_f, time, steps):
        Efetua a operação de mutação em um membro da população.
    '''
    def __init__(self, thetas_i, thetas_f, time, steps):
        '''
        Parâmetros
        ----------
        theta_i:
            Lista de ângulos inicias das juntas.
        theta_i:
            Lista de ângulos finais das juntas.
        time: float
            Tempo de duração da simulação.
        steps: int
            Lista com número de pontos presentes para divisão 
            e criação de novos polinômios.
        
        '''
        times, thetas, omegas = [], [], []
        angles_elements = [list(create_angles(theta_i, theta_f, time, steps)) for theta_i, theta_f in zip(thetas_i, thetas_f)]
        for i in range(5):
            deltas, delta_thetas, delta_omegas = create_angles(thetas_i[i], thetas_f[i], time, steps)
            times.append(deltas)
            thetas.append(delta_thetas)
            omegas.append(delta_omegas)
        self.joint = [[times[i], thetas[i], omegas[i]] for i in range(5)]

    def find_points(self, theta_i, time, steps):    
        '''
        Método para criação das curvas de cada elemento, aplicação de
        cinemática para encontrar variáveis cartesianas, verificação de 
        métrica e constraints
        
        Parâmetros
        ----------
        theta_i: list
            Lista de ângulos inicias das juntas.
        time: float
            Tempo de duração da simulação.
        steps: int
            Lista com número de pontos presentes para divisão 
            e criação de novos polinômios.
        '''
        mentor = Mentor()
        points, polynomies = [], []
        for j in range(steps-1):
            sub_polynomies = []
            for i in range(5):
                sub_polynomies.append(Polinomy(self.joint[i][0][j], self.joint[i][0][j+1], self.joint[i][1][j], self.joint[i][1][j+1], self.joint[i][2][j], self.joint[i][2][j+1], number=POINTS))
            polynomies.append(sub_polynomies)
        angles = [np.array([theta_i[i]]) for i in range(5)]
        for i in range(steps-1):
            for j in range(len(angles)):
                angles[j] = np.concatenate((angles[j], polynomies[i][j].thetas[1:])) 
        index = np.min([angles[0].shape[0], angles[1].shape[0], angles[2].shape[0], angles[3].shape[0], angles[4].shape[0]])
        self.angle = [angles[i][0:index] for i in range(5)]
        angles = np.transpose(self.angle)
        dist = 0 
        for i, angle in enumerate(angles):
            new_pos, rot = mentor.get_position(angle, z_axis=5)
            if i == 0:
                old_pos = new_pos
            else:
                dist+=np.sqrt((old_pos[0]-new_pos[0])**2+(old_pos[1]-new_pos[1])**2+(old_pos[2]-new_pos[2])**2)
                old_pos = new_pos
            points.append(new_pos[0:3])
        self.dist = dist
        self.points = points
        self.constraint = False
        self.test_velocity(time)

    def test_velocity(self, time):
        '''
        Método para verificação dos constraints com base na velocidade
        máxima permitida definida nas constantes.

        Parâmetros
        ----------
        time: float
            Tempo de duração da simulação.
        '''
        for i in range(self.angle[0].shape[0]-1):
            for j in range(5):
                if abs(self.angle[j][i+1]-self.angle[j][i])/(time/self.angle[j].shape[0]) > MAXIMUM_VELOCITY:
                    self.constraint = True
                    break
            if self.constraint:
                break                
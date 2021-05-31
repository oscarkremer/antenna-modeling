'''
Este script executa os cálculos de cinemática direta e inversa do 
robô mentor. Primeiramente inseri-se um conjunto de ângulos de juntas
os quais são utilizados na cinemática direta para cálculo da posição
e orientação. Com as variáveis do atuador final aplica-se cinemática 
inversa, permitindo verificar se os ângulos encontrados são iguais aos
inseridos.

Com o ambiente mentor ativado no conda este script pode ser 
executado com o comando make kinematics.
'''
import numpy as np
from src.mentor import Mentor
from src.utils.input import input_angles

if __name__ == "__main__":
    robot = Mentor()
    angles  = input_angles()
    pos, rot = robot.get_position(angles)
    print('Vetor de Posição: ')
    print(pos)
    print('Matriz de Rotação: ')
    print(rot)
    print('Thetas: {}'.format(180*np.array(robot.get_angles(pos,rot))/np.pi))
    print('Confirmalção Cinemática Inversa: {}'.format(robot.get_position(robot.get_angles(pos,rot))))
'''
Este script define a classe Mentor, incluindo métodos
matemáticos para cálculo de cinématica direta e inversa
do manipulador, orientação e posição do atuador final,
verificação de posição e correções de quadrante. 

Este arquivo pode ser importado como um módulo utilizando:
from src.mentor import Mentor
'''
import itertools
import numpy as np

class Mentor:
    '''
    Classe de um robo mentor. A classe será definida
    por elementos que obedeçam as leis de cinemática do robô mentor,
    além de limitações (constraints) de velocidade e posição.

    Atributos
    ----------
    alpha : list
        List contendo tempo, angulos e velocidades discretizadas.
    a: list
        Métrica de distância percorrida do elemento.
    d: list
        Coordenadas cartesianos do elemento ao longo do tempo.

    Métodos
    -------
    test_inverse_kinematics(self, pos, rot, tag_theta1=True, tag_theta2=True, tag_theta3=True):  
        Método de teste de todos os ângulos, para verificação se os quadrantes se encontram corretos.
    get_angles(self, pos, rot):
        Cálculo de cinemática inversa para todas as possibilidades de quadrantes entre os três primeiros
        atuadores rotacionais do robô.
    _inverse_kinematics(self, pos, orientation, tag_theta1=True, tag_theta2=True, tag_theta3=True):
        Cálculo da cinemática inversa para uma cenário específico.
    verify(self, pos, rot, returned_pos, returned_rot):
        Verificação se posição encontrada pela cinemática é compatível com o desejado.
    get_position(self, theta, z_axis=0):
        Aplicação de cinemática direta para encontrar posição do atuador.
    get_orientation(self, alpha, beta, gamma):
        Aplicação de cinemática direta para encontrar orientação do atuador.    
    separate(self, matrix, z_axis):
        Separação da matriz extendida nas parcelas de rotação e translação.    
    denavit(self, theta, i):
        Aplicação de denavit hartenberg para encontrar matriz de transformação.
    fix_quadrante(self, sin, cos):
        Ajustar quadrante a partir de sin e cos de ângulo.
    fix_theta(self, param, tag, angle):
        Ajustar quadrante a partir de qualquer parâmetro de entrada possivel.
    '''
    def __init__(self):
        '''
        Não há parâmetros de entrada neste método. Apenas a definição dos
        parâmetros de Denavit-Hartenberg.
        '''
        self.alpha = [0, -np.pi/2, 0, 0, -np.pi/2]
        self.a = [0, 0, 17.2739, 15.5, 0]
        self.d = [0, 0, 0, 0, 0]

    def test_inverse_kinematics(self, pos, rot, tag_theta1=True, tag_theta2=True, tag_theta3=True):  
        '''
        Método de teste de todos os ângulos, fazendo o cálculo da cinemática inversa e 
        a verificação se os quadrantes estão corretos. Ajustes dos ângulos com base nas 
        tags especificadas. A função de cinemática inversa é chamada para encontrar os ângulos, 
        para então serem passados para função verify que testará se os mesmos estão corretos.

        Parâmetros
        ----------
        pos: list
            Lista com informações de posição do atuador final do manipuldor 
            no espaço cartesiano.
        rot: list
            Matriz de rotação do atuador final do manipulador no espaço 
            cartesiano.
        tag_theta1: bool(opcional):
            Tag de verificação se há erro no quadrante do ângulo 1.
        tag_theta2: bool(opcional):
            Tag de verificação se há erro no quadrante do ângulo 2.
        tag_theta3: bool(opcional):
            Tag de verificação se há erro no quadrante do ângulo 3.

        Retorna
        -------
        tag:
            Verificação de há erro ou não na cinemática inversa.
        theta:
            Lista com ângulos calculados.
        '''

        theta = self._inverse_kinematics(pos, rot, tag_theta1=tag_theta1, tag_theta2=tag_theta2, tag_theta3=tag_theta3)
        returned_pos, returned_rot = self.get_position(theta)        
        tag = self.verify(pos, rot, returned_pos, returned_rot)
        return tag, theta

    def get_angles(self, pos, rot):
        '''
        Método para efetuar cinemática inversa em todos os cenários possíveis
        de erros de quadrantes. Com base nisso é possível identificar se a 
        posição/orientação inserida pelo usuário é atingível ou não pelo manipulador.

        Parâmetros
        ----------
        pos: list
            Lista com informações de posição do atuador final do manipuldor 
            no espaço cartesiano.
        rot: list
            Matriz de rotação do atuador final do manipulador no espaço 
            cartesiano.

        Retorna
        -------
        tag:
            Booleano que identifica (se verdadeiro) se há erro e a posição é impossível
        theta:
            Lista com ângulos calculados para cinemática inversa.
        '''
        for element in itertools.product([True, False],[False, True],[False, True]):
            tag, theta = self.test_inverse_kinematics(pos, rot, element[0], element[1], element[2])
            if tag:
                return False, theta
        return True, theta 

    def _inverse_kinematics(self, pos, orientation, tag_theta1=True, tag_theta2=True, tag_theta3=True):
        '''
        Método para cálculo da cinemática inversa. Esta função implementa o equacionamento
        mostrado na documentação do projeto. A correção de quadrante é feita com base
        nas tags de cada ângulo.

        Parâmetros
        ----------
        pos: list
            Lista com informações de posição do atuador final do manipuldor 
            no espaço cartesiano.
        rot: list
            Matriz de rotação do atuador final do manipulador no espaço 
            cartesiano.
        tag_theta1: bool(opcional):
            Tag de verificação se há erro no quadrante do ângulo 1.
        tag_theta2: bool(opcional):
            Tag de verificação se há erro no quadrante do ângulo 2.
        tag_theta3: bool(opcional):
            Tag de verificação se há erro no quadrante do ângulo 3.
        
        Retorna
        -------
        theta:
            Lista com ângulos calculados para cinemática inversa.
        '''
        theta = []
        theta1 = np.nan_to_num(self.fix_theta(pos, tag_theta1, 'theta1'))
        theta.append(theta1)
        theta3 = self.fix_theta((pos[0]**2+pos[1]**2+pos[2]**2-self.a[2]**2-self.a[3]**2)/(2*self.a[2]*self.a[3]), tag_theta3, 'theta3')
        theta3 = np.nan_to_num(theta3)
        theta2  = np.nan_to_num(-theta3+self.fix_theta(((-np.cos(theta3)*self.a[2]-self.a[3])*pos[2]+np.sin(theta3)*self.a[2]*(np.sin(theta1)*pos[1]+pos[0]*np.cos(theta1)))/(np.power(np.cos(theta1)*pos[0]+np.sin(theta1)*pos[1], 2)+pos[2]*pos[2]), tag_theta2, 'theta2'))
        theta.append(theta2)
        theta.append(theta3)
        sin_theta4 = np.sin(theta2+theta3)*orientation[2][2] - np.cos(theta1)*np.cos(theta2+theta3)*orientation[0][2] - np.sin(theta1)*np.cos(theta2+theta3)*orientation[1][2]
        cos_theta4 = -np.cos(theta2+theta3)*orientation[2][2] - np.cos(theta1)*np.sin(theta2+theta3)*orientation[0][2] - np.sin(theta1)*np.sin(theta2+theta3)*orientation[1][2]
        theta.append(self.fix_quadrante(sin_theta4, cos_theta4))
        sin_theta5 = np.sin(theta1)*orientation[0][0] - np.cos(theta1)*orientation[1][0]
        cos_theta5 = np.sin(theta1)*orientation[0][1] - np.cos(theta1)*orientation[1][1]
        theta.append(self.fix_quadrante(sin_theta5, cos_theta5))
        return theta

    def verify(self, pos, rot, returned_pos, returned_rot):
        '''
        Método para verificação se os ângulos encontrados realmente
        atingem a posição solicitada.

        Parâmetros
        ----------
        pos: list
            Lista com informações de posição correta do atuador final do manipuldor 
            no espaço cartesiano.
        rot: list
            Matriz de rotação correta do atuador final do manipulador no espaço 
            cartesiano.
        returned_pos: list
            Lista com informações de posição do atuador final do manipuldor 
            no espaço cartesiano.
        returned_rot: list
            Matriz de rotação do atuador final do manipulador no espaço 
            cartesiano.
        
        Retorna
        -------
        Retorna boolean para identifição se a posição está ou não correta.
        '''
        diff = pos[:3]-returned_pos[:3]
        rot_norm = np.linalg.norm(rot-returned_rot)
        if rot_norm > 0.001 or abs(diff[0])>0.001 or abs(diff[1]>0.001) or abs(diff[2]>0.001):
            return False
        else:
            return True

    def get_position(self, theta, z_axis=0):
        '''
        Método para encontrar posição a partir da cinemática direta.

        Parâmetros
        ----------
        theta: Numpy Array
            Vetor com os ângulos das juntas.
        z_axis: float (opcional)
            Distância do sistema de coordenadas do atuador final em 
            relação à última junta rotacional do Mentor.

        Retorna
        -------
        Matrizes de rotação e orientação que definem o sistema de coordenada
        final do atuador no sistema cartesiano da base.
        '''
        matrix = np.matmul(self.denavit(theta, 3),self.denavit(theta, 4))
        for i in range(3):
            matrix = np.matmul(self.denavit(theta, 2-i), matrix)
        return self.separate(matrix, z_axis)
    
    def get_orientation(self, alpha, beta, gamma):
        '''
        Método para verificação dos constraints. Com base na velocidade
        máxima permitida definida nas constantes.

        Parâmetros
        ----------
        alpha: float
            Ângulo de orientação do eixo X.
        beta: float
            Ângulo de orientação do eixo Y.
        gamma: float
            Ângulo de orientação do eixo Z.            .

        Retorna
        -------
        Matriz de rotação encontrada a partir dos ângulos alpha-beta-gamma.
        '''
        return [[np.cos(alpha)*np.cos(beta), np.cos(alpha)*np.sin(beta)*np.sin(gamma) - np.sin(alpha)*np.cos(gamma), np.cos(alpha)*np.sin(beta)*np.cos(gamma) + np.sin(alpha)*np.sin(gamma) ],
                        [np.sin(alpha)*np.cos(beta), np.sin(alpha)*np.sin(beta)*np.sin(gamma) + np.cos(alpha)*np.cos(gamma), np.sin(alpha)*np.sin(beta)*np.cos(gamma) - np.cos(alpha)*np.sin(gamma) ],
                        [-np.sin(beta), -np.cos(beta)*np.sin(gamma), np.cos(beta)*np.cos(gamma)]]

    def separate(self, matrix, z_axis):
        '''
        Método para separação da matriz de Denavit-Hartenberg 
        da última junta nas parcelas de posição e orientação.

        Parâmetros
        ----------
        theta: Numpy Array
            Vetor com os ângulos das juntas.
        z_axis: float (opcional)
            Distância do sistema de coordenadas do atuador final em 
            relação à última junta rotacional do Mentor.

        Retorna
        -------
        pos: float
            Matriz 3x3 de rotação do robô que define a orientação da garra.
        np.array(rot):
            Matriz 3x3 de rotação do robô que define a orientação da garra.
        '''
        pos =  np.matmul(matrix, [0,0,z_axis,1])
        rot = [matrix[i][:-1] for i in range(3)]
        return pos, np.array(rot)

    def denavit(self, theta, i):
        '''
        Método para criação de matriz de Denavit-Hartenberg de uma
        determinada transição de sistemas coordenados.

        Parâmetros
        ----------
        theta: list
            Lista com ângulos das juntas em determinado instante.
        i: int
            Index que definirá entre quais sistemas coordenados está acontecendo 
            está transformação de coordenadas.

        Retorna
        -------
        Matrix 4x4 extendida e triangular por blocos que contém 
            matriz de rotação (3x3) e vetor de posição (3x1).
        '''
        return  [[np.cos(theta[i]), -np.sin(theta[i]),  0,   self.a[i]],
            [np.sin(theta[i])*np.cos(self.alpha[i]), np.cos(theta[i])*np.cos(self.alpha[i]), -np.sin(self.alpha[i]), -np.sin(self.alpha[i])*self.d[i]],
            [np.sin(theta[i])*np.sin(self.alpha[i]), np.cos(theta[i])*np.sin(self.alpha[i]), np.cos(self.alpha[i]), np.cos(self.alpha[i])*self.d[i]],
            [0, 0, 0, 1]]

    def fix_quadrante(self, sin, cos):
        '''
        Método para correção de quadrante a partir do par sin/cos.

        Parâmetros
        ----------
        sin: float
            Tempo de duração da simulação.
        cos: float
            Tempo de duração da simulação.

        Retorna
        -------
        Ângulo ajustado para o quadrante correto
        '''
        if sin >= 0 and cos >= 0:
            return np.arcsin(abs(sin))
        if sin >= 0  and cos <= 0:
            return np.pi - np.arcsin(abs(sin))
        if sin <= 0 and cos >= 0:
            return 2*np.pi - np.arcsin(abs(sin))
        if sin <= 0 and cos <= 0:
            return np.arccos(abs(cos)) + np.pi

    def fix_theta(self, param, tag, angle):
        '''
        Método correção dos três primeiros ângulos do robô.

        Parâmetros
        ----------
        parm: float
            Parâmetro que será analisado para o ajuste, podendo ser, 
            lista de posição (x,y), sin e cos.
        tag: boolean
            Tag de verificação que determina se o ângulo está correto ou não.
        angle: str
            Descritivo sobre qual ângulo é, permitindo vincular com 
            operação matemática especifica.

        Retorna
        -------
        Ângulo corrigido para o quadrante correto.
        '''
        if angle=='theta1':
            if abs(param[0]) < 0.001:  
                param = np.pi/2
            else:
                param = np.arctan(param[1]/ param[0])
        if tag and param>=0:
            if angle == 'theta1':
                return param
            if angle == 'theta2':
                return np.arcsin(abs(param))       
            if angle == 'theta3':
                return np.arccos(abs(param))       
        if not tag and param>=0:
            if angle == 'theta1':
                return np.pi + param
            if angle == 'theta2':
                return np.pi-np.arcsin(abs(param))
            if angle == 'theta3':
                return 2*np.pi-np.arccos(abs(param))   
        if tag and param<=0:
            if angle == 'theta1':
                return 2*np.pi - abs(param)
            if angle == 'theta2':
                return 2*np.pi - np.arcsin(abs(param))
            if angle == 'theta3':   
                return np.pi-np.arccos(abs(param))       
        if not tag and param<=0:
            if angle == 'theta1':
                return np.pi - abs(param)
            if angle == 'theta2':
                return np.pi + np.arcsin(abs(param))
            if angle == 'theta3':
                return np.pi + np.arccos(abs(param))
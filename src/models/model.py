'''
Este script define a classe da população do algoritmo genético.

Este arquivo pode ser importado como um módulo utilizando:
from src.models import Population
'''
import itertools
import random
import numpy as np
import pandas as pd
from tqdm import tqdm
from .node.node import Node
from src.utils.numerical import create_angles
from src.utils.constants import PATH

class Population:
    '''
    Classe da população de robôs mentor's que são otimizados
    quanto à distância percorrida em uma determinada trajetória.

    Atributos
    ----------
    size : int
        Inteiro que representa o tamanho da população ao longo
        da execução do processo de otimização.
    initial: int
        Inteiro que representa o tamanho inicial da população.
    generations: int
        Número de iterações que serão feitas.
    p_c: float
        Probabilidade de ocorrer o processo de cross-over dentro 
        de um novo elemento da população.
    p_m: float
        Probabilidade de ocorrer o processo de mutação dentro de 
        um novo elemento da população.
    mode: srt
        String que define o modo de cálculo das probabilidades de 
        cross-over e mutação, caso seja adaptativo as probabilidades 
        serão calculadas com base no comportamento médio da população.

    Métodos
    -------
    initialization(self, theta_i, theta_f, time, steps):
        Inicializa a população e cria os primeiros polinômios de 
        acordo com o número de passos definidos.

    generation(self, population,  theta_i, theta_f, time, steps):
        Efetua as operações de cross-over e mutação em cima de uma 
        determinada geração.

    save_csv(self, population, actual_best, best_of_all):
        Salva as informações de trajetória e desempenhos em um .csv

    prob_adaptation(self, fitness):
        Efetua a adaptação das probabilidades de acordo com o desempenho médio

    analysis(self, population):
        Calcula os parâmetros estatísticos presentes e utilizados na adaptação.

    selection(self, population, number_bests = 2):
        Efetua a seleção natural nos elementos da geração.

    cross_over(self, parents, theta_i, theta_f, time, steps):
        Efetua a operação de cross-over com base em dois pais de um membro.

    mutation(self, member, theta_i, theta_f, time, steps):
        Efetua a operação de mutação em um membro da população.
    '''

    def __init__(self, initial, size, generations, p_c, p_m, mode='adaptive'):
        '''
        Parâmetros
        ----------
        size: int
            Tamanho da população ao longo das gerações
        initial: int
            Tamanho da população na primeira geração
        generations: int
            Número de gerações que serão simuladas
        p_c: float
            Probabilidade de cross-over, para o adaptivo este valor representa o 
            fator constante da função de adaptação
        p_m: float
            Probabilidade de mutação, para o adaptivo este valor representa o 
            fator constante da função de adaptação
        mode: str (opcional)
            Modo de tratamento das probabilidades de mutação e cross-over, 
            default é 'adaptive'
        '''
        self.size = size
        self.initial = initial
        self.generations = generations
        self.p_c = p_c
        self.p_m = p_m
        self.mode = mode

    def initialization(self, theta_i, theta_f, time, steps):
        '''
        Método para inicialização da primeira geração, 
        criando o número de elementos definidos com o atributo
        initial. Com o constraint dos elementos verificado
        considera-se apenas os membros válidos.

        Parâmetros
        ----------
        theta_i: list
            Lista de ângulos inicias da junta.
        theta_f: list
            Lista de ângulos finais que serão ocupados pelo robô.
        time: list
            Lista de instantes de tempo.
        steps: list
            Lista com número de pontos presentes para divisão 
            e criação de novos polinômios.

        Retorna
        -------
        population: list
            Lista de elementos e métrica, que contém cada elemento
            com seu respectivo desempenho.
        '''
        population = []
        i = 0
        while i < self.initial:
            element = Node(theta_i, theta_f, time, steps)
            element.find_points(theta_i, time, steps)
            if not element.constraint:
                i+=1
                metric = element.dist
                population.append([metric, element])
        return population
    
    def generation(self, population,  theta_i, theta_f, time, steps):
        '''
        Método para operação dos manipuladores genéticos em uma deter-
        minada população. As operação são executadas dentro de um loop
        para n iterações, onde n é o número de gerações. A melhor traje-
        tória e o melhor comportamento da população ao longo das gerações
        são salvos em um arquivo .csv

        Parâmetros
        ----------
        population: list
            Lista de elements e métrica, que contém cada elemento
            com seu respectivo desempenho.
        theta_i: list
            Lista de ângulos inicias da junta.
        theta_f: list
            Lista de ângulos finais que serão ocupados pelo robô.
        time: list
            Lista de com o vetor de instantes de tempo.
        steps: list
            Lista com número de pontos presentes para divisão 
            e criação de novos polinômios.
        '''
        best_of_generation = []
        actual_best = []
        for i in tqdm(range(self.generations)):
            population = self.selection(population, number_bests = self.size)
            actual_best.append(self.selection(population, number_bests = 1)[0][0])
            self.analysis(population)
            members = [member[1] for member in population]
            combinations = list(itertools.product(members, repeat=2))
            for combination in combinations:
                new_element = self.cross_over(combination)
                new_element.find_points(theta_i, time, steps)
                if not new_element.constraint:
                    population.insert(len(population), [new_element.dist, new_element]) 
            self.analysis(population)
            for member in population:
                mutation = self.mutation(member, theta_i, theta_f, time, steps)
                mutation.find_points(theta_i, time, steps)
                if not mutation.constraint:
                    population.append([mutation.dist, mutation])
            best_of_generation.append(self.selection(population, number_bests = 1)[0][0])
            best_generation = self.selection(population, number_bests = 1)[0]
        best_of_all = self.selection(population, number_bests = 1)[0]
        self.save_csv(actual_best, best_of_all)


    def save_csv(self, actual_best, best_of_all):
        '''
        Método para passar as informações do melhor elemento encontrado
        e do desempenho da população ao longo do processo de otimização
        para dataframes no formato .csv.

        Parâmetros
        ----------
        actual_best: list
            Lista contendo melhor desempenho da população ao longo
            das gerações, salvando a informação em points.csv.
        best_of_all: list
            Lista contendo a métrica e o melhor elemento, utilizada 
            para salvar a trajetória no espaço cartesiano e de juntas
            no arquivo points.csv.
        '''
        dataframe, points = pd.DataFrame(), pd.DataFrame()
        dataframe['distance'] = actual_best
        x, y, z = [], [], []
        for point in best_of_all[1].points:
            x.append(point[0])
            y.append(point[1])
            z.append(point[2])          
        points['x'] = x
        points['y'] = y
        points['z'] = z
        for i in range(5):
            points['theta{}'.format(i+1)] = best_of_all[1].angle[i]
        points.to_csv('data/results/points.csv', index=False)
        dataframe.to_csv('data/results/results.csv', index=False)

    def prob_adaptation(self, fitness):
        '''
        Método para adaptação das probabilidades de mutação e cross-over, 
        a otimização é feita com base no desempenho máximo e médio da 
        população.

        Parâmetros
        ----------
        fitness: float
            Variável que contém (distancia percorrida)^-1, buscando a minimização
            desta métrica.
        '''
        if self.mode == 'adaptive':
            if fitness > self.average and (self.maximum-self.average) > 0.00001:
                self.p_c = (self.maximum-fitness)/(self.maximum-self.average)
                self.p_m = (self.maximum-fitness)/(self.maximum-self.average)
            else:
                self.p_c = 0.7
                self.p_m = 0.8

    def analysis(self, population):
        '''
        Método para cálculo dos parâmetros estatísticos vinculados
        ao processo de adaptação das probabilidades, calcula-se 
        desvio padrão, valor máximo e médio.

        Parâmetros
        ----------
        population: list
            Lista de elementos e métrica, que contém cada elemento
            com seu respectivo desempenho.
        '''
        fitness = [1/member[0] for member in population]
        self.std = np.std(np.array(fitness))
        self.maximum = max(fitness)
        self.average = sum(fitness)/len(fitness)

    def selection(self, population, number_bests = 2):
        '''
        Método para seleção dos melhores elementos de um conjunto 
        de membros, população. 

        Parâmetros
        ----------
        population: list
            Lista de elements e métrica, que contém cada elemento
            com seu respectivo desempenho.
        number_bests: integer (opcional)
            Número de melhores elementos a serem selecionados 
            (default=2)

        Retorna
        -------
        Lista dos melhores elementos devidamente ordenados, em ordem 
        decrescente pela métrica de fit.
        '''
        return (sorted(population, key = getitem))[0:number_bests]
   
    def cross_over(self, parents):
        '''
        Método para operação de cross-over em dois elementos pais, 
        criando um novo elemento para futura geração.

        Parâmetros
        ----------
        parents: list
            Lista contendo dois elementos que formarão um novo membro 
            para a próxima geração.
        
        Retorna
        -------
        Novo elemento criado pela troca do material genético presente
        entre os dois pais. Lista de elements e métrica que contém cada elemento
        com seu respectivo desempenho.
        
        '''
        self.prob_adaptation(1/parents[0].dist)
        for i in range(5):
            if random.random() <= self.p_c:
                parents[0].joint[i] = parents[1].joint[i]
        return parents[0]

    def mutation(self, member, theta_i, theta_f, time, steps):
        '''
        Método para operação de mutação dentro de um elemento.

        Parâmetros
        ----------
        member: list
            Lista contendo um elemento e sua métrica.
        theta_i: list
            Lista de ângulos inicias da junta.
        theta_f: list
            Lista de ângulos finais que serão ocupados pelo robô.
        time: list
            Lista de com o vetor de instantes de tempo
        steps: list
            Lista com número de pontos presentes para divisão 
            e criação de novos polinômios.

        Retorna
        -------
        Novo elemento fruto da mutação.
        '''
        element = Node(theta_i, theta_f, time, steps)
        for i in range(5):
            element.joint[i] = member[1].joint[i]
        self.prob_adaptation(1/member[0])
        for i in range(5):
            if random.random() <= self.p_m:
                deltas, delta_thetas, delta_omegas = create_angles(theta_i[i], theta_f[i], time, steps)
                member[1].joint[i] = [deltas, delta_thetas, delta_omegas]
        return member[1]

def getitem(item):
    return item[0]
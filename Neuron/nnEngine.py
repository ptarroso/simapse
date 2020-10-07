#!/usr/bin/env python
'''
SIMAPSE - simulation maps for ecological niche modelling
Version 2.00
Copyright (C) 2010  Pedro Tarroso

Please cite: 
"Tarroso, P., Carvalho, S. & Brito, J.C. (2012) Simapse - Simulation
Maps for Ecological Niche Modelling. Methods in Ecology and Evolution
doi: 10.1111/j.2041-210X.2012.00210.x"

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

from random import random
from math import exp

class NN():
    def __init__(self, scheme=[], iterations=1000, LR=0.9, momentum=0.0, verbosity = 1):
        # l - representa layer
        # n - representa neuronio
        # w - representa weight
        
        self.scheme = scheme
        self.structure(scheme) #Create network

        self.iterations = iterations
        self.LearningRate = float(LR)
        self.momentum = float(momentum)

        self.Patterns = None #O LoadData calcula as Patterns
        self.trainInputs = None
        self.trainOutputs = None
        self.currentPat = None
        self.errPat = None
        self.GlobalError = None

        self.func = sigm
        self.dfunc = dsigm

        self.verbosity = verbosity


    def trainnet(self, verbose = None):
        '''Trains the network with loaded inputs.
           Verbose level:
                0 - Nothing is printed
                1 - Prints Iteration number | Network error'''

        if verbose == None:
            verbose = self.verbosity

        for i in range(0, self.iterations):
            for p in range(0, self.Patterns):
                self.currentPat = p

                #calculate the current network output
                #and error for this pattern
                #FEED FORWARD
                self.feedforward()
                # change weights on network
                # BACK PROPAGATION
                self.backpropag()

            #error for this iteration
            if verbose == 1:
                self.neterror(self.trainInputs, self.trainOutputs, 'SSerror')
                print("iteration = {} | RMS error = {}".format(i, self.GlobalError))

    def loaddata(self, inputs, targets):
        #Le os dados para o treino
        self.trainInputs = inputs
        self.trainOutputs = targets
        self.Patterns = len(self.trainInputs)
        #self.rndWeights() # Randomization of weights must be manual
        self.__somestats()

    def __somestats(self):
        '''Calculate some statistics about inputs'''
        inputs = self.trainInputs
        nvar = self.scheme[0]
        n = len(inputs)

        #variable average
        average = [0.0] * nvar
        for i in range(n):
            average = [average[x] + inputs[i][x] for x in range(nvar)]
        average = [average[x] / n for x in range(nvar)]

        #variable standard deviation
        stdev = [0.0] * nvar
        for i in range(n):
            stdev = [stdev[x] + (inputs[i][x] - average[x])**2 for x in range(nvar)]
        stdev = [stdev[x] / n for x in range(nvar)]

        #Assign self variables
        self.average = average
        self.stdev = stdev   

    def testnet(self, inputs, verbose = None):
        '''Tests the network with a sequence of inputs and returns
           the predictided output.
           verbose level:
                0 - nothing is printed
                1 - prints \'pattern number | original target | network output\''''

        if verbose == None:
            verbose = self.verbosity

        #Accepts a single list of inputs or a list of sequence of inputs
        if type(inputs[0]) is not list:
            inputs = [inputs]
        self.trainInputs = inputs
        nInputs = self.Patterns = len(inputs)
        finalresults = []
        for p in range(nInputs):
            self.currentPat = p
            self.feedforward()
            result = self.values[-1][:]
            finalresults.append(result) 
            if verbose == 1:
                print("Pattern = {} | {} | predicted = {}".format(p+1,self.trainOutputs[p], result))

        return finalresults

    def neterror(self, inputs = None, targets = None, errorType = 'SSerror'):
        '''Calculates the overall error of the network.
           Options for error type are:
           RMSerror - Root Mean Square error (default)
           SSerror  - Sum of Squared error'''
        try:
            if errorType not in ['RMSerror', 'SSerror']:
                raise SyntaxError('Error type must be \'RMSerror\' or \'SSerror\'!')

            if inputs is not None and targets is not None:
                self.loaddata(inputs, targets)

            if errorType == 'RMSerror':
                self.__RMSerror()
            elif errorType == 'SSerror':
                self.__SSerror()
        
            error = self.GlobalError
            return error

        except (SyntaxError) as e:
            print(e)

    def __RMSerror(self):
        '''Calculates the Root Mean Square Error os the network.'''
        nOutputs = self.scheme[-1]
        temp = [0.0] * nOutputs
        for p in range(self.Patterns):
            self.currentPat = p #NAO ESQUECER QUE E APENAS UM INDICE DA LISTA
            self.feedforward()
            temp = [temp[x] + (self.errPat[x])**2 for x in range(nOutputs)]
        RMSerror = [(temp[x] / self.Patterns)**0.5 for x in range(nOutputs)]
        self.GlobalError = RMSerror

    def __SSerror(self):
        '''Calculates the Sum of Squared Errors of the network.'''
        nOutputs = self.scheme[-1]
        temp = [0.0] * nOutputs
        for p in range(self.Patterns):
            self.currentPat = p #NAO ESQUECER QUE E APENAS UM INDICE DA LISTA
            self.feedforward()
            temp = [temp[x] + (self.errPat[x])**2 for x in range(nOutputs)]
        SSerror = [0.5 * (temp[x]) for x in range(nOutputs)]
        self.GlobalError = SSerror        

    def backpropag(self):
        nlayers = len(self.scheme) - 1 #Without input layer
        LR = self.LearningRate
        M = self.momentum
        curPat = self.currentPat
        scheme = self.scheme
        dfunc = self.dfunc
        nlayers = len(scheme)
        weights = self.weights
        values = self.values
        changes = self.changes
        tInputs = self.trainInputs[curPat][:]
        out_errors = self.errPat[:]           # Erros da layer seguinte - comeca pelos erros de output

        #Para o output
        errors = [0.0] * (scheme[-2] + 1)     #Erros para a ultima hidden layer com BIAS
        for n in range(scheme[-1]):           #Para cada neuronio de output
            change = 0.0
            derivative = dfunc(values[-1][n])
            delta = derivative * out_errors[n]        #TODO o feedforward armazena a derivativa em self.derivatives! Posso utilizar aqui?
            for w in range(scheme[-2]):       #Para cada weight associado ao neuronio (SEM BIAS)
                errors[w] += weights[-1][n][w] * out_errors[n]  #Calculo do erro para a ultima hidden layer
                change = LR * delta * values[-2][w] + M * changes[-1][n][w]
                weights[-1][n][w] = weights[-1][n][w] - change
                changes[-1][n][w] = change
            errors[-1] += weights[-1][n][-1] * out_errors[n]
            change =  LR * delta                        #Calculo do BIAS
            weights[-1][n][-1] = weights[-1][n][-1] - change
            
        #Para as hidden layers
        for l in range(nlayers - 2, 1, -1):                     #Para todas as hidden layers (sem output e sem a priemira com ligacao aos inputs)
            prevL_errors = [0.0] * (scheme[l-1] + 1)        #erros para a previous layer (next layer in the reversed sequence) + 1 para BIAS
            for n in range(scheme[l]):                    #Para todos os neuronios da layer l
                change = 0.0
                derivative = dfunc(values[l-1][n]) # -1 porque a values nao tem inputs
                delta = derivative * errors[n] 
                for w in range(scheme[l-1]):         #Para todos os weights do neuronio n da layer l SEM BIAS (igual ao numero de neuronios da layer anterior)
                    prevL_errors[w] += weights[l-1][n][w] * errors[n]
                    change = LR * delta * values[l-2][w] + M * changes[l-1][n][w]
                    weights[l-1][n][w] = weights[l-1][n][w] - change
                    changes[l-1][n][w] = change
                prevL_errors[-1] += weights[l-1][n][-1] * errors[n]
                change = LR * delta
                weights[l-1][n][-1] = weights[l-1][n][-1] - change
            errors = prevL_errors[:]                            #Copia os erros para a seuqencia seguinte

        # Para a primeira hidden layer
        for n in range(scheme[1]):                             #Neuronios da primeira hidden layer
            change = 0.0
            derivative = dfunc(values[0][n])
            delta = derivative * errors[n]
            for w in range(scheme[0]):                         #Numero de inputs (sem o bias)
                change = LR * delta * tInputs[w] + M * changes[0][n][w]
                weights[0][n][w] = weights[0][n][w] - change
                changes[0][n][w] = change
            change = LR * delta
            weights[0][n][-1] = weights[0][n][-1] - change


    def feedforward(self):
        #Mudou de neuralnet() para feedforward()
        #calculate the outputs of the hidden and output neurons
        #the hidden neurons are tanh
        scheme = self.scheme
        hValues = self.values
        weights = self.weights
        curPat = self.currentPat
        nlayers = len(hValues) #POSSO PASSAR PARA O __INIT__? 
        tInputs = self.trainInputs[curPat][:]
        tInputs.append(1.0) #BIAS
        derivatives = self.derivatives
        func = self.func
        dfunc = self.dfunc

        #FEEDFORWARD
        # Para os inputs - 1a hidden layer
        for n in range(scheme[1]):      #Para cada neuronio da primeira layer
            hValues[0][n] = 0.0         #Reset aos valores para o somatorio
            for w in range(scheme[0]): #Para cada weight da primeira layer sem o BIAS
                hValues[0][n] += tInputs[w] * weights[0][n][w]
            hValues[0][n] += weights[0][n][-1]          # Soma o weight do BIAS
            hValues[0][n] = func(hValues[0][n])
            derivatives[0][n] = dfunc(hValues[0][n])
        
        # Para as hidden layers
        for l in range(1, nlayers):      #Para todas as hidden layers e output
            for n in range(scheme[l+1]): #Para todos os neuronios da hLayer l (+1 devido ao input)
                hValues[l][n] = 0.0     #Resete do valor do neuronio para o somatorio
                for w in range(scheme[l]): #Para todos os weights menos BIAS
                    hValues[l][n] += hValues[l-1][w] * weights[l][n][w]
                hValues[l][n] += weights[l][n][-1]
                hValues[l][n] = func(hValues[l][n])
                derivatives[l][n] = dfunc(hValues[l][n])
        self.calcError() # Tirei do Testnet, assim sempre q faz 1 feedforward calcula o erro
        self.values = hValues

    def calcError(self):
        #calculate the error for the outputs
        #TODO Analisar o codigo!!!
        #TODO Pode passar para o backpropagation??? Assim posso usar a funcao para testar 
        #sem fazer calculos desnecessarios
        #Mas deve ser necessario para calcular o overallerror
        curPat = self.currentPat
        hValues = self.values
        nOutputs = self.scheme[-1]
        errPat = [hValues[-1][o] - self.trainOutputs[curPat][o]for o in range(nOutputs)]
        # nao deve ser preciso a atribuicao depois
        self.errPat = errPat

    def pderiv(self, inputs):
        '''Processes the partial derivatives of the output vs. input
           Inputs must be [[i1],[i2],[i3],...], or [inputs]'''
        if len(inputs) >= 1 and type(inputs[0]) is list:
            pderiv = []
            for i in inputs:
                pderiv.append(self.__pderiv(i))
        elif len(inputs) == self.scheme[0] and type(inputs[0]) is float:
            pderiv = self.__pderiv(inputs)
        elif len(inputs) < 1 or type(inputs[0]) is not list:
            print("Inputs must be [[i1],[i2],[i3],...] or [inputs]")
        return pderiv

    def _weight_gen(self, layer, rep = 1, remove_bias = True):
        '''Returns the weights of the layer sorted by the
           sequence of the neurons in the previous layer
           neurons. It repeats the sequence 'rep' number
           of times and it optionally removes the BIAS weight'''
        bias = int(bool(remove_bias))
        weights = self.weights
        new_weights = []
        w = weights[layer]
        for r in range(rep):
            for n_previous in range(len(w[0]) - bias):
                for n_next in range(len(w)):
                    new_weights.append(w[n_next][n_previous])
        return new_weights

    def __pderiv(self, inputs):
        '''Returns the partial derivatives of the network output
           with respect to the input. This function only processes
           one sequence of inputs by sprawling the network.'''
        self.trainInputs = [inputs]
        self.currentPat = 0
        self.feedforward()

        deriv = self.derivatives
        weights = self.weights
        scheme = self.scheme[1:] # Without input
        nlayers = len(scheme)
        noutputs = scheme[-1]

        pderiv = [[0.0 for y in range(len(inputs))] for x in range(len(deriv[-1]))]
        for i in range(len(inputs)):
            # variable product has all the weights connected to input node i
            product = [weights[0][x][i] for x in range(len(weights[0]))]

            for l in range(nlayers-1): # output layer removed
                # n is the number of connections available to between current layers
                # m is number of nodes with derivatives in layer l
                n, m = len(product), len(deriv[l])
                # each deriv node value for layer l is repeated 
                # by the number of connections it has (the same as the number of 
                # neurons of previous layer)
                new_deriv = deriv[l] * int(n / m)
                # the products are multiplied by the derivatives sequentially
                product = [product[x] * new_deriv[x] for x in range(n)]
                # w has all the weights from next layer sorted by the neurons of
                # the previous layer and repeated for the number of connections
                # of the current calculation
                w = self._weight_gen(l+1, n)
                # new product is created by repeating each value by the number of
                # neuron in the next layer
                sprawl_product = repeat(product, len(deriv[l+1]))
                # sprawled product is multiplied by the sequence of weights between
                # current and next layer
                product = [sprawl_product[x] * w[x] for x in range(len(sprawl_product))]

            n, m = len(product), len(deriv[-1])
            new_deriv = deriv[-1] * int(n / m)
            product = [product[x] * new_deriv[x] for x in range(n)]

            for o in range(noutputs):
                for x in range(0, n, noutputs):
                    pderiv[o][i] += product[x+o]
        return pderiv

    def derivexp(self, inputs):  
        '''May be used for simple networks with several inputs, one
           hidden layer and one output. The algorithm follows the
           simplified formula of Dimopoulos, 1995.
           Just used to confirm the partial derivatives result from
           more complex network given by __pderiv() 
           MAY BE DELETED!!!'''

        self.trainInputs = [inputs]
        self.currentPat = 0
        self.feedforward()

        deriv = self.derivatives
        weights = self.weights
        scheme = self.scheme[1:] # Without input
        nlayers = len(scheme)
        noutputs = scheme[-1]
        pderiv = []
        for i in range(len(inputs)):
            value = 0
            for ni in range(scheme[0]):
                #weights e layer 0 que tem weights entre inputs e 1a hidden
                #e layer 1 que tem entre 1a hidden e output
                #pderiv e a layer 0, que corresponde aos  neuronios
                #da 1a hidden layer
                value += weights[0][ni][i] * deriv[0][ni] * weights[1][0][ni]
            value = value * deriv[-1][0]
            pderiv.append(value)
        return pderiv
            

    def structure(self, network):
        '''Creates the structure of the network
           Weights are representes as the weights that connect
           to a given neuron. A BIAS is added to each weight 
           sequence for a neuron.'''
        try:
            if self.scheme == []:
                raise SyntaxError('Network scheme cannot be empty!')

            nlayers = len(network)
            values, derivatives, weights, changes = [], [], [], []

            for layer in range(1, nlayers):
                values.append([])
                derivatives.append([])
                weights.append ([])
                changes.append([])
                nneurons, p_nneurons = network[layer], network[layer-1]+1 # +1 para o BIAS nos weights
                i = layer-1
                for neuron in range(nneurons):
                    values[i].append(0.0)
                    derivatives[i].append(0.0)
                    weights[i].append([0.0]* (p_nneurons ))
                    changes[i].append([0.0]* (p_nneurons))          

            self.values = values
            self.derivatives = derivatives
            self.weights = weights
            self.changes = changes

        except (SyntaxError) as e:
            print(e)

    def rndWeights(self):
        #Inicializa todos os weights com numeros aleatorios entre [-0.5, 0.5]
        w = self.weights
        for l in range(len(w)):
            for n in range(len(w[l])):
                w[l][n] = [(random() - 0.5) for x in w[l][n]]
        self.weights = w

    def XORexample(self):
        #Fazer um data loader com calculo automatico do self.pattern
        # o __init__ corre logo esta funcao
        #TODO O numero de patterns tem de ser igual ao numero de outputs
        print("Loading XOR example...\n")
        Inputs = [[1, 0], [0, 1], [1, 1], [0, 0]]
        Targets = [[1.0], [1.0], [0.0], [0.0]]
        self.loaddata(Inputs, Targets)
        self.rndWeights()
        self.trainnet(verbose = 1)

        for p in range(0, self.Patterns):
            self.currentPat = p #NAO ESQUECER QUE E APENAS UM INDICE DA LISTA
            self.feedforward()
            print("Pattern = {} | real = {} | predicted = {}".format(p+1, self.trainOutputs[p], self.values[-1]))
        derivs = self.pderiv(Inputs)
        for i in range(len(Inputs)):
            print('Input: {} | Partial derivative: {}'.format(Inputs[i], derivs[i]))

def repeat(lst, rep):
    '''Creates a new list where each item is repeated 'rep'
       number of times, mantaining the original order:
       repeat([1,2,3], 2) = [1,1,2,2,3,3]'''
    new_lst = []
    for item in lst:
        for r in range(rep):
            new_lst.append(item)
    return new_lst

def savenet(net, netfile):
    '''Saves the trained network to a file.'''
    f = open(netfile, 'w')
    netvars = vars(net)
    for var in netvars:
        value = netvars[var]
        if hasattr(value, '__call__'):
            value = value.__name__
        line = str(var) + ';' + str(value) + '\n' 
        f.write(line)
    f.close()
    
def loadnet(netfile):
    '''Loads a saved trained network.'''
    f = open(netfile, 'r')
    data = f.readlines()
    f.close()
    netvars = {}
    for line in data:
        line = line.split(';')
        netvars[line[0]] = eval(line[1])
    net = NN(netvars['scheme'])
    net.__dict__ = netvars
    return net            
        
def tanh(x):
    result = None        
    if x > 20:
        result = 1
    elif x < -20:
        result = -1
    else:
        a = exp(x)
        b = exp(x * -1)
        result = (a-b) / (a+b)
    return result

def dtanh(y):
    result = 1 - (y**2)
    return result

def sigm(x):
    e = 2.7182818284590451
    if x < -700: #avoid overflow on large networks
        return 1 / (1 + e**700)
    else:
        return 1 / (1 + e**(-x))

def dsigm(y):
    return y * (1 - y)


if __name__ == '__main__':
    nn = NN([2,3,1], iterations=10000, LR=0.8, momentum=0.0)
    nn.XORexample()

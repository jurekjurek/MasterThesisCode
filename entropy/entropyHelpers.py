'''
Function to go through the test diagrams one by one 
and modify them accordingly. 
'''
import numpy as np
from qiskit.circuit import Parameter
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.converters import dag_to_circuit, circuit_to_dag

def SetUpDiags(circuit: QuantumCircuit, parameterDict: dict, WordsToForget: list = ['woman', 'person', 'meal'], oneQubit = True, show: bool = True):
    '''
    Given a circuit as an Input, we remove the noun parameters belonging to the first noun in the sentence. 

    We then, depending on the shape of the circuit, replace the noun parameters by an entangled Bell state or we forget about it. 

    In the end, we replace the verb parameters by the parameters learned by the model. 

    The reason we are doing this is to introduce a model with real uncertainty. 

    wordsToForget: list of strings of words that we want to forget ('man', 'woman', 'person', ...)

    '''
    # remove noun parameters and set the rest of the parameters 

    # start by extracting circuit data 
    circuitData = circuit.data

    if show:     
        print('circuit that we start with: ')
        print(circuit)


    # list to store the nounParameters to be removed 
    indexList = []

    # boolean that tells us if the noun of interest is daggered or not 
    dagger = False 

    for i in range(len(circuitData)): 

        # check if this gate has a parameter of a word we want to forget about
        for word in WordsToForget:
            if word in str(circuitData[i].operation.params):
                # store the parameter in a list
                indexList.append(i)

                print('we forget the following word: ', word)

                qubitOfInterest = circuitData[i].qubits[0].index

                if '†' in str(circuitData[i].operation.params):
                    dagger = True 
                    
                # break


    # FIRST: UPDATE THE PARAMETERS, THEN REMOVE THE GATES FROM THE CIRCUIT 

    # get parameters of the remaining circuit 
    params = circuit.parameters


    # set the parameters of the remaining words based on the models weights 
    for i in range(len(params)):
        # print(str(params[i]))
        # print(parameterDict[str(params[i])])
        print(params[i])
        if str(params[i]) in parameterDict:
            circuit = circuit.assign_parameters({params[i]: (parameterDict[str(params[i])]/(2*np.pi))})
        else: 
            print('TESTSETESTSETS')
            # if the model only learned parameters for the word that are not daggered
            if '†' in str(params[i]):
                if oneQubit == True: 
                    tempParam = str(params[i])
                    # reverse string 
                    tempParam = tempParam[::-1]

                    # remove elements 2,3,4,5
                    tempParam = tempParam[:2] + 'n__' + tempParam[5:]

                    tempParam = tempParam[::-1]
                    print(tempParam)             

                    circuit = circuit.assign_parameters({params[i]: -(parameterDict[tempParam]/(2*np.pi))}) 
                if oneQubit == False:
                    tempParam = str(params[i])
                    # reverse string 
                    tempParam = tempParam[::-1]

                    # remove elements 2,3,4,5
                    tempParam = tempParam[:2] + 'n__' + tempParam[5:]

                    tempParam = tempParam[::-1]
                    print(tempParam)             

                    circuit = circuit.assign_parameters({params[i]: (parameterDict[tempParam]/(2*np.pi))}) 

    if show: 
        print('Circuit wiht new parameters: ')
        print(circuit)
    

    # remove the parameters of the words we want to forget
    circuit.data.pop(indexList[0])
    circuit.data.pop(indexList[1]-1)
    circuit.data.pop(indexList[2]-2)

    if show: 
        print('circuit with removed words to forget: ')
        print(circuit)


    # if the word we want to forget about is an effect, we simply forget about the qubit 
    if dagger: 
        # remove measurement 
        circuitData = circuit.data

        for i in range(len(circuitData)):

            if str(circuitData[i].operation.name) == 'measure': 
                if circuitData[i].qubits[0].index == qubitOfInterest: 

                    # remove measurement
                    circuit.data.pop(i)
                    break
        
        # that's it 
        

    # if the word we want to forget about is a state, we introduce a bell state and clip it to the beginning of the circuit 
    if not dagger: 
        # introduce circuit with 5 qubits to be composed with the initial circuit 
        bellCircuit = QuantumCircuit(circuit.num_qubits + 1)
        # bell state on first two qubits 
        bellCircuit.h(0)

        # entangle the artificial qubit with the qubit of interest
        bellCircuit.cx(0,qubitOfInterest + 1)

        # now, add one 'artificial' qubit to the initial circuit 
        q = QuantumRegister(1, 'q')
        circuit = circuit.reverse_bits()
        circuit.add_bits(q)
        circuit = circuit.reverse_bits()

        # now that the two circuits are equal in length, combine them 
        circuit = bellCircuit.compose(circuit) #, qubits=[0,1,2,3,4])

        # now, we have our circuit 

    if show: 
        print('circuit with discarded qubit: ')
        print(circuit)

    circuit.draw(output= 'mpl', filename='testtest.png')


    return circuit, not dagger


from lambeq import BobcatParser
from lambeq import RemoveCupsRewriter
from lambeq import IQPAnsatz

from pytket.extensions.qiskit import tk_to_qiskit





'''
Now, for the density matrices 
'''
from qiskit import QuantumCircuit
from qiskit.quantum_info import DensityMatrix, partial_trace
from qiskit.visualization import plot_state_city



def GetDensityMatrix(circuit: QuantumCircuit, traceList: list, addedQubit: bool):
    '''
    Given a QuantumCircuit representing a sentence, this function returns the density matrix of the *sentence* qubit 
    ''' 
    import qiskit.quantum_info as qi 
    from qiskit.providers.aer import AerSimulator
    
    dmSimulation = AerSimulator()

    job = dmSimulation.run(circuit.decompose(reps = 1))



    result = job.result()

    # print(result)

    data = result.data()

    rho = data.get('density_matrix')

    # print('CIRCUIT IN DENSITY MATRIX')
    # print(circuit)

    # differentiate between dagger and not dagger (one has one qubit less)
    # dagger
    # if circuit.num_qubits == 4: 
    #     sentenceDM = partial_trace(state=rho, qargs=[0, 2, 3])

    # # not dagger 
    # elif circuit.num_qubits == 5: 
    #     sentenceDM = partial_trace(state=rho, qargs=[0, 1, 2, 4])

    # # simple sentence like man prepares food, without any adjectives, man will be daggered there 
    # elif circuit.num_qubits == 3: 
    #     sentenceDM = partial_trace(state=rho, qargs=[0,2])

    if not addedQubit:
        sentenceDM = partial_trace(state=rho, qargs= traceList)
    if addedQubit: 
        traceList = [i+1 for i in traceList]
        traceList.insert(0,0)
        # print('modified tracelist if appended:', traceList)
        sentenceDM = partial_trace(state=rho, qargs= traceList)

    # else: 
    #     print('number of qubits: ', circuit.num_qubits)
    #     print(circuit)

    return sentenceDM

'''
now for the application of the code
'''

def Main(listOfCircuits: list, parameterDict: dict, wordsToForget: list):
    '''
    Iterate over all circuits in the test_circuits list and for each of them create the density matrix as explained above 
    ''' 

    listOfDensityMatrices: list = []


    from pytket.extensions.qiskit import tk_to_qiskit

    numQubitsList: list = []

    for i in range(len(listOfCircuits)): 

        print('iteration number: ', i)

        tempCirc = listOfCircuits[i]

        # if i == 8: 
        #     tempCirc.draw()

        tempCircTK = tempCirc.to_tk()

        tempCircQiskit = tk_to_qiskit(tempCircTK)

        # if i == 8: 
        #     print(tempCircQiskit)

        circuitData = tempCircQiskit.data

        # list of qubits that shall be forgotten when traced over
        traceList = []
        for i in range(len(circuitData)):
            if str(circuitData[i].operation.name) == 'measure':
                traceList.append(circuitData[i].qubits[0].index)
        # print('traceList: ', traceList)

        alteredTempCirc, addedQubit = SetUpDiags(tempCircQiskit, parameterDict, wordsToForget)

        # to access density matrix later on 
        alteredTempCirc.save_density_matrix() 

        densityMatrix = GetDensityMatrix(alteredTempCirc, traceList, addedQubit) 
        
        listOfDensityMatrices.append(densityMatrix)

        numQubitsList.append(tempCircQiskit.num_qubits)

    return listOfDensityMatrices, numQubitsList
    




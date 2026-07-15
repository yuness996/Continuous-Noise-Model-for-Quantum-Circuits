import numpy as np

from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, pauli_error

L0 = ['00000','10010','01001','10100','01010','11011','00110','11000','11101','00011','11110','01111','10001','01100','10111','00101'] # Superposed states of |0L>
L1 = ['11111','01101','10110','01011','10101','00100','11001','00111','00010','11100','00001','10000','01110','10011','01000','11010'] # Superposed states of |1L>
def post_processing_decoding(dic):
    """After measuring all of the qubits we post-process the result.
    This function will interpret whether the measured state is in which logic state (|0L> or |1L>) for 1 or 2 encoded qubits.
    It returns a dictionary in |0>/|1> basis."""
    if len(next(iter(dic))) < 10:# 1 encoded qubit with & w/o error correction
        dic_out = {'0':0,'1':0}
        for item, count in dic.items():
            iteminv = item[::-1]
            if iteminv[:5][::-1] in L0:
                dic_out['0'] += count
            elif iteminv[:5][::-1] in L1:
                dic_out['1'] += count
    return dic_out

def runcircuitnoisy(qc,err,p_err,shots):
    """This function runs a quantum circuit qc in a Pauli noise channel.
    err : array of the type of Pauli errors (X,Y,Z).
    p_err : array with the probability of each error.
    shots : number of shots."""
    arr = []
    for i in range(len(err)):
        arr.append((err[i], p_err[i]))
    arr.append(('I', 1 - np.sum(p_err)))
    error = pauli_error(arr)
    noise_model = NoiseModel()
    noise_model.add_all_qubit_quantum_error(error, ["u1", "u2", "u3"])
    sim_noise = AerSimulator(noise_model=noise_model)
    qc_noise = transpile(qc, sim_noise)
    result = sim_noise.run(qc_noise,shots=shots).result()
    counts = result.get_counts(0)
    return counts

def runcircxtimes(qc,err,p_err,shots,times):
    """Function that runs qc with shots number of shots on a noisy simulator and repeats the simulation 
    times number of times then returns the the mean value of the searched value and the error bar."""
    prob = []
    for i in range(int(times)):
        counts = runcircuitnoisy(qc,err,p_err,int(shots))
        countsL = post_processing_decoding(counts)
        prob.append(countsL)
    return prob

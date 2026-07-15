import numpy as np

from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit

from Oracles import *
from PostProcessing import *


def single_encoded_qubit(p_err, shots, times, correction, n_h):
    if correction:
        fn = 6
    else:
        fn = 0

    qq = QuantumRegister(7 + fn, "q")
    cc = ClassicalRegister(7 + fn, "c")

    qc = QuantumCircuit(qq, cc)

    qc.append(encode(), range(7))

    for _ in range(n_h):
        qc.barrier()
        qc.append(Logical_Hadamard(), range(7))
        qc.barrier()

    if correction:
        qc.append(detection(), range(13))

        qc = decode_correction(
            qc,
            qq,
            [cc[7], cc[8], cc[9], cc[10], cc[11], cc[12]],
        )

    for i in range(7):
        qc.measure(qq[i], cc[i])

    return runcircxtimes(qc,["X", "Y", "Z"],[p_err] * 3,shots,times)


def Run_AllNoise(a, b, num, shots, times, correction, n_h):
    """
    Run the [[7,1,3]] encoded circuit over a logarithmic noise range.
    """
    p_errors = np.logspace(a, b, num=num)

    arr = []

    for pe in p_errors:
        arr.append(single_encoded_qubit(p_err=pe / 3,shots=shots,times=times,correction=correction,n_h=n_h))

    return arr


def save_table(table, filename):
    """
    Save the simulation output to a text file.

    Each row corresponds to one physical error rate.
    Each row contains `times` decoded logical-count dictionaries.
    """
    with open(filename, "a", encoding="utf-8") as file:
        for row in table:
            file.write(" ".join(map(str, row)) + "\n")


if __name__ == "__main__":
    save_table(
        Run_AllNoise(
            a=-6,
            b=-1,
            num=20,
            shots=1000,
            times=10,
            correction=True,
            n_h=10,
        ),
        "test.txt",
    )

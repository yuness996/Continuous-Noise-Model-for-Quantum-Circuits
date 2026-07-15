from qiskit import *
import numpy as np

from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, transpile
from qiskit.circuit import QuantumCircuit, Qubit, Clbit


def encode():
    """
    Encoding circuit for the [[7,1,3]] Steane code.

    The circuit prepares one logical qubit encoded into
    seven physical qubits.
    """

    qc = QuantumCircuit(7, name="Encoding")

    # Prepare the three ancilla superposition states.
    qc.h([0, 1, 2])

    # Apply the CNOT network defining the Steane encoding.
    qc.cx(6, 4)
    qc.cx(6, 5)

    qc.cx(0, 3)
    qc.cx(0, 5)
    qc.cx(0, 6)

    qc.cx(1, 3)
    qc.cx(1, 4)
    qc.cx(1, 6)

    qc.cx(2, 3)
    qc.cx(2, 4)
    qc.cx(2, 5)

    return qc


def detection():
    """
    Syndrome-extraction circuit for the [[7,1,3]] code.

    Qubits 0-6 are data qubits.
    Qubits 7-12 are syndrome ancillas.

    The first three ancillas measure X stabilizers.
    The last three ancillas measure Z stabilizers.
    """

    qc = QuantumCircuit(13, name="Detection")

    # Prepare syndrome ancillas in |+>.
    qc.h(range(7, 13))

    # X stabilizer measurements.
    qc.cx(7, 0)
    qc.cx(7, 2)
    qc.cx(7, 4)
    qc.cx(7, 6)

    qc.cx(8, 0)
    qc.cx(8, 1)
    qc.cx(8, 4)
    qc.cx(8, 5)

    qc.cx(9, 0)
    qc.cx(9, 1)
    qc.cx(9, 2)
    qc.cx(9, 3)

    # Z stabilizer measurements.
    qc.cz(10, 0)
    qc.cz(10, 2)
    qc.cz(10, 4)
    qc.cz(10, 6)

    qc.cz(11, 0)
    qc.cz(11, 1)
    qc.cz(11, 4)
    qc.cz(11, 5)

    qc.cz(12, 0)
    qc.cz(12, 1)
    qc.cz(12, 2)
    qc.cz(12, 3)

    # Rotate ancillas back to the computational basis.
    qc.h(range(7, 13))

    return qc


def Logical_Hadamard():
    """
    Logical Hadamard gate.

    For the Steane code the logical Hadamard is transversal,
    i.e. it consists of a Hadamard applied to every data qubit.
    """

    qc = QuantumCircuit(7, name="$H_L$")

    qc.barrier()
    qc.h(range(7))
    qc.barrier()

    return qc

def decode_correction(qc,qbits,cbits):
    """
    Measure the syndrome ancillas and apply error correction for the
    [[7,1,3]] Steane code.

    The six measured syndrome bits are interpreted using a hard-decision
    lookup table implemented with nested conditional statements. Each
    syndrome is mapped to the corresponding single-qubit Pauli correction.

    Parameters
    ----------
    qc : QuantumCircuit
        Circuit containing the syndrome extraction stage.

    qbits : QuantumRegister
        Quantum register containing the seven data qubits followed by
        the six syndrome ancillas.

    cbits : list[Clbit]
        Classical bits used to store the measured syndrome.

    Returns
    -------
    QuantumCircuit
        Circuit with syndrome measurement, conditional corrections,
        and ancilla reset operations.
    """
    for i in range(6):
        qc.measure(qbits[i+7],cbits[i])
    with qc.if_test((cbits[0], 0)) as else_0:#0---
        with qc.if_test((cbits[1], 0)) as else_1:#00--
            with qc.if_test((cbits[2], 0)) as else_2:#000-
                with qc.if_test((cbits[3], 0)) as else_3:#0000-
                    with qc.if_test((cbits[4], 0)) as else_4:#00000-
                        with qc.if_test((cbits[5], 1)):#000001
                            qc.x(qbits[3])
                    with else_4:#00001-
                        with qc.if_test((cbits[5], 0)) as else_5:#000010
                            qc.x(qbits[5])
                        with else_5:#000011
                            qc.x(qbits[1])
                with else_3:#0001-
                    with qc.if_test((cbits[4], 0)) as else_6:#00010-
                        with qc.if_test((cbits[5], 0)) as else_7:#000100
                            qc.x(qbits[6])
                        with else_7:#000101
                            qc.x(qbits[2])
                    with else_6:#00011-
                        with qc.if_test((cbits[5], 0)) as else_8:#000110
                            qc.x(qbits[4])
                        with else_8:#000111
                            qc.x(qbits[0])
            with else_2:#001-
                with qc.if_test((cbits[3], 0)):#0010-
                    with qc.if_test((cbits[4], 0)):#00100-
                        with qc.if_test((cbits[5], 0)) as else_9:#001000
                            qc.z(qbits[3])
                        with else_9:#001001
                            qc.y(qbits[3])
        with else_1:#01--
            with qc.if_test((cbits[2], 0)) as else_10:#010-
                with qc.if_test((cbits[3], 0)):#0100-
                    with qc.if_test((cbits[4], 0)) as else_11:#01000-
                        with qc.if_test((cbits[5], 0)):#010000
                            qc.z(qbits[5])
                    with else_11:#01001-
                        with qc.if_test((cbits[5], 0)):#010010
                            qc.y(qbits[5])
            with else_10:#011-
                with qc.if_test((cbits[3], 0)):#0110-
                    with qc.if_test((cbits[4], 0)) as else_12:#01100-
                        with qc.if_test((cbits[5], 0)):#011000
                            qc.z(qbits[1])
                    with else_12:#01101-
                        with qc.if_test((cbits[5], 1)):#011011
                            qc.y(qbits[1])
    with else_0:#1---
        with qc.if_test((cbits[1], 0)) as else_13:#10--
            with qc.if_test((cbits[2], 0)) as else_14:#100-
                with qc.if_test((cbits[3], 0)) as else_15:#1000-
                    with qc.if_test((cbits[4], 0)):#10000-
                        with qc.if_test((cbits[5], 0)):#100000
                            qc.z(qbits[6])
                with else_15:#1001-
                    with qc.if_test((cbits[4], 0)):#10010-
                        with qc.if_test((cbits[5], 0)):#100100
                            qc.y(qbits[6])
            with else_14:#101-
                with qc.if_test((cbits[3], 0)) as else_16:#1010-
                    with qc.if_test((cbits[4], 0)):#10100-
                        with qc.if_test((cbits[5], 0)):#101000
                            qc.z(qbits[2])
                with else_16:#1011-
                    with qc.if_test((cbits[4], 0)):#10110-
                        with qc.if_test((cbits[5], 1)):#101101
                            qc.y(qbits[2])
        with else_13:#11--
            with qc.if_test((cbits[2], 0)) as else_17:#110-
                with qc.if_test((cbits[3], 0)) as else_18:#1100-
                    with qc.if_test((cbits[4], 0)):#11000-
                        with qc.if_test((cbits[5], 0)):#110000
                            qc.z(qbits[4])
                with else_18:#1101-
                    with qc.if_test((cbits[4], 1)):#11011-
                        with qc.if_test((cbits[5], 0)):#110110
                            qc.y(qbits[4])
            with else_17:#111-
                with qc.if_test((cbits[3], 0)) as else_19:#1110-
                    with qc.if_test((cbits[4], 0)):#11100-
                        with qc.if_test((cbits[5], 0)):#111000
                            qc.z(qbits[0])
                with else_19:#1111-
                    with qc.if_test((cbits[4], 1)):#11111-
                        with qc.if_test((cbits[5], 1)):#111111
                            qc.y(qbits[0])
    for i in range(6):
        qc.reset(qbits[i+7])
    return qc

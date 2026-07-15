"""Noisy circuit blocks for the [[7,1,3]] Steane code."""

import numpy as np
from qiskit import QuantumCircuit

# ---------------------------------------------------------------------
# Continuous coherent-error model
# ---------------------------------------------------------------------

def r_error_gate(qc, theta, phi, qubit):
    """Apply the sampled single-qubit coherent error rotation."""
    qc.z(qubit)
    qc.rx(2*theta, qubit)
    qc.rz(-2*phi, qubit)
    qc.z(qubit)
    return qc
def insert_error(qc, pos, k):
    """
    Insert independent Gaussian coherent errors on the selected qubits.

    Each angle is sampled from N(0, (k/3)^2).
    """
    for p in pos:
        theta = np.random.normal(loc=0.0, scale=k/3)
        phi = np.random.normal(loc=0.0, scale=k/3)
        qc = r_error_gate(qc, theta, phi, p)
    return qc
# ---------------------------------------------------------------------
# [[7,1,3]] encoding circuit
# ---------------------------------------------------------------------

def encode(k):
    """Construct the noisy encoding circuit for one Steane-code block."""
    qc = QuantumCircuit(7,name='Encoding')
    qc.h([0,1,2])
    qc = insert_error(qc,[0,1,2],k)
    qc.cx(6,4)
    qc.cx(6,5)
    qc.barrier()
    qc.cx(0,3)
    qc.cx(0,5)
    qc.cx(0,6)
    qc.barrier()
    qc.cx(1,3)
    qc.cx(1,4)
    qc.cx(1,6)
    qc.barrier()
    qc.cx(2,3)
    qc.cx(2,4)
    qc.cx(2,5)
    qc.barrier()
    return qc
# ---------------------------------------------------------------------
# Syndrome-extraction circuit
# ---------------------------------------------------------------------

def detection(k):
    """
    Construct the noisy six-ancilla syndrome-extraction circuit.

    Data qubits are 0--6 and syndrome ancillas are 7--12.
    Noise is inserted after each single-qubit Hadamard layer.
    """
    qc = QuantumCircuit(13,name='Detection')
    qc.h(range(7,13))
    qc = insert_error(qc,range(7,13),k)
    qc.barrier()
    qc.cx(7,0)
    qc.cx(7,2)
    qc.cx(7,4)
    qc.cx(7,6)
    qc.barrier()
    qc.cx(8,0)
    qc.cx(8,1)
    qc.cx(8,4)
    qc.cx(8,5)
    qc.barrier()
    qc.cx(9,0)
    qc.cx(9,1)
    qc.cx(9,2)
    qc.cx(9,3)
    qc.barrier()
    qc.h([0,2,4,6])
    qc = insert_error(qc,[0,2,4,6],k)
    qc.cx(10,0)
    qc.cx(10,2)
    qc.cx(10,4)
    qc.cx(10,6)
    qc.h([0,2,4,6])
    qc = insert_error(qc,[0,2,4,6],k)
    qc.barrier()
    qc.h([0,1,4,5])
    qc = insert_error(qc,[0,1,4,5],k)
    qc.cx(11,0)
    qc.cx(11,1)
    qc.cx(11,4)
    qc.cx(11,5)
    qc.h([0,1,4,5])
    qc = insert_error(qc,[0,1,4,5],k)
    qc.barrier()
    qc.h([0,1,2,3])
    qc = insert_error(qc,[0,1,2,3],k)
    qc.cx(12,0)
    qc.cx(12,1)
    qc.cx(12,2)
    qc.cx(12,3)
    qc.h([0,1,2,3])
    qc = insert_error(qc,[0,1,2,3],k)
    qc.barrier()
    qc.h(range(7,13))
    qc = insert_error(qc,range(7,13),k)
    return qc
# ---------------------------------------------------------------------
# Logical two-block CNOT
# ---------------------------------------------------------------------

def cnot():
    """Return the transversal logical CNOT between two Steane blocks."""
    qc = QuantumCircuit(14,name='CNOT')
    for i in range(7):
        qc.cx(i,i+7)
    return qc
# ---------------------------------------------------------------------
# Logical Hadamard
# ---------------------------------------------------------------------

def h_L(k):
    """Construct the noisy transversal logical Hadamard circuit."""
    qc = QuantumCircuit(7,name='$H_L$')
    qc.h(range(7))
    qc = insert_error(qc,[0,1,2,3],k)
    return qc
# ---------------------------------------------------------------------
# Syndrome decoding and recovery
# ---------------------------------------------------------------------

def decode_correction(qc, qbits, cbits):
    """
    Measure the six syndrome ancillas and apply hard-decision recovery.

    The nested conditional tree maps each six-bit syndrome to the
    corresponding single-qubit X, Y, or Z correction.
    """
    # Measure the six syndrome ancillas.
    for i in range(6):
        qc.measure(qbits[i + 7], cbits[i])
    
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
    
    # Reset the ancillas so that they can be reused.
    for i in range(6):
        qc.reset(qbits[i + 7])
    return qc
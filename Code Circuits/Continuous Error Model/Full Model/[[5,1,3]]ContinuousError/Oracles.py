"""Continuous coherent-noise circuits for the [[5,1,3]] code."""

import time

import matplotlib.pyplot as plt
import numpy as np

from qiskit import (
    ClassicalRegister,
    QuantumCircuit,
    QuantumRegister,
    transpile,
)
from qiskit.circuit import Clbit, Qubit
from qiskit.circuit.library import UGate
from qiskit.visualization import (
    circuit_drawer,
    latex,
    plot_bloch_multivector,
    plot_histogram,
)
from qiskit_aer import AerSimulator
from qiskit_aer.noise import (
    NoiseModel,
    QuantumError,
    ReadoutError,
    depolarizing_error,
    pauli_error,
    thermal_relaxation_error,
)


# ---------------------------------------------------------------------
# Continuous coherent-error model
# ---------------------------------------------------------------------

def r_error_gate(qc, theta, phi, qubit):
    """
    Apply the continuous single-qubit error rotation.

    The error is implemented through the gate sequence

        Z R_x(2 theta) R_z(-2 phi) Z.

    Parameters
    ----------
    qc : QuantumCircuit
        Circuit to which the error is added.

    theta : float
        Random rotation-angle error.

    phi : float
        Random phase-angle error.

    qubit : int or Qubit
        Target qubit.

    Returns
    -------
    QuantumCircuit
        Circuit containing the error rotation.
    """

    qc.z(qubit)
    qc.rx(2 * theta, qubit)
    qc.rz(-2 * phi, qubit)
    qc.z(qubit)

    return qc


def insert_error(qc, pos, k):
    """
    Insert independent coherent errors on a set of qubits.

    For each qubit in ``pos``, two independent angles are sampled from
    zero-mean Gaussian distributions with standard deviation ``k / 3``.

    Parameters
    ----------
    qc : QuantumCircuit
        Circuit to modify.

    pos : iterable
        Qubit indices on which errors are inserted.

    k : float
        Noise-strength parameter. The Gaussian standard deviation is
        taken as k / 3.

    Returns
    -------
    QuantumCircuit
        Circuit containing the sampled coherent errors.
    """

    for p in pos:

        # Alternative von Mises model used in earlier tests:
        #
        # theta = np.random.vonmises(mu=0.0, kappa=9 / k**2)
        # phi = np.random.vonmises(
        #     mu=0.0,
        #     kappa=9 / k**2,
        # ) % (2 * np.pi)

        # Sample independent small-angle Gaussian errors.
        theta = np.random.normal(
            loc=0.0,
            scale=k / 3,
        )

        phi = np.random.normal(
            loc=0.0,
            scale=k / 3,
        )

        # Add the sampled unitary error to the selected qubit.
        qc = r_error_gate(
            qc,
            theta,
            phi,
            p,
        )

    return qc


# ---------------------------------------------------------------------
# [[5,1,3]] encoding circuit
# ---------------------------------------------------------------------

def encode(k):
    """
    Construct the noisy encoding circuit for the [[5,1,3]] code.

    Coherent errors are inserted after each layer of single-qubit gates.

    Parameters
    ----------
    k : float
        Continuous-noise strength.

    Returns
    -------
    QuantumCircuit
        Five-qubit encoding circuit.
    """

    qc = QuantumCircuit(
        5,
        name="$Encoding$",
    )

    # First single-qubit layer.
    qc.h([4, 3, 2, 1])
    qc.z(0)

    # Add noise after the single-qubit gates.
    qc = insert_error(
        qc,
        range(5),
        k,
    )

    # First entangling layer.
    qc.cx(4, 0)
    qc.cx(3, 0)
    qc.cx(2, 0)
    qc.cx(1, 0)

    # Second single-qubit layer.
    qc.h([0, 3, 1])

    qc = insert_error(
        qc,
        [0, 3, 1],
        k,
    )

    # Second entangling layer.
    qc.cx(4, 0)
    qc.cx(4, 3)
    qc.cx(2, 1)

    # Third single-qubit layer.
    qc.h([3, 1, 2, 0])

    qc = insert_error(
        qc,
        [3, 1, 2, 0],
        k,
    )

    # Additional Hadamard layer on qubit 0.
    qc.h(0)

    qc = insert_error(
        qc,
        [0],
        k,
    )

    # Final entangling layer.
    qc.cx(3, 2)
    qc.cx(1, 0)

    # Final single-qubit layer.
    qc.h([2, 0])

    qc = insert_error(
        qc,
        [2, 0],
        k,
    )

    return qc


# ---------------------------------------------------------------------
# Logical Hadamard circuit
# ---------------------------------------------------------------------

def Logical_Hadamard(k):
    """
    Construct the noisy logical Hadamard circuit for the [[5,1,3]] code.

    The logical operation consists of transversal physical Hadamard gates
    followed by a qubit permutation implemented using SWAP gates.

    Parameters
    ----------
    k : float
        Continuous-noise strength.

    Returns
    -------
    QuantumCircuit
        Logical Hadamard circuit acting on five physical qubits.
    """

    qc = QuantumCircuit(
        5,
        name="$H_L$",
    )

    # Apply a physical Hadamard to each data qubit.
    qc.h(range(5))

    # Add coherent noise after the Hadamard layer.
    qc = insert_error(
        qc,
        range(5),
        k,
    )

    # Apply the qubit permutation required by the logical Hadamard.
    qc.swap(0, 1)
    qc.swap(1, 3)
    qc.swap(1, 4)

    return qc


# ---------------------------------------------------------------------
# Error-detection circuit
# ---------------------------------------------------------------------

def detection(k):
    """
    Construct the noisy syndrome-extraction circuit for the [[5,1,3]] code.

    The circuit uses four syndrome ancillas. Each ancilla measures one
    stabilizer generator through a five-qubit subcircuit.

    Data qubits:
        0, 1, 2, 3, 4

    Syndrome ancillas:
        5, 6, 7, 8

    Coherent errors are inserted only after single-qubit Hadamard gates.

    Parameters
    ----------
    k : float
        Continuous-noise strength.

    Returns
    -------
    QuantumCircuit
        Nine-qubit syndrome-extraction circuit.
    """

    # --------------------------------------------------------------
    # First stabilizer-measurement subcircuit
    # --------------------------------------------------------------

    cs1 = QuantumCircuit(5)

    qc_h_qubits = [1, 2]

    cs1.h(qc_h_qubits)
    cs1 = insert_error(
        cs1,
        qc_h_qubits,
        k,
    )

    cs1.cx(4, 0)
    cs1.cx(4, 1)
    cs1.cx(4, 2)
    cs1.cx(4, 3)

    cs1.h(qc_h_qubits)
    cs1 = insert_error(
        cs1,
        qc_h_qubits,
        k,
    )

    # --------------------------------------------------------------
    # Second stabilizer-measurement subcircuit
    # --------------------------------------------------------------

    cs2 = QuantumCircuit(5)

    qc_h_qubits = [1, 2]

    cs2.h(qc_h_qubits)
    cs2 = insert_error(
        cs2,
        qc_h_qubits,
        k,
    )

    cs2.cx(4, 0)
    cs2.cx(4, 1)
    cs2.cx(4, 2)
    cs2.cx(4, 3)

    cs2.h(qc_h_qubits)
    cs2 = insert_error(
        cs2,
        qc_h_qubits,
        k,
    )

    # --------------------------------------------------------------
    # Third stabilizer-measurement subcircuit
    # --------------------------------------------------------------

    cs3 = QuantumCircuit(5)

    qc_h_qubits = [0, 1]

    cs3.h(qc_h_qubits)
    cs3 = insert_error(
        cs3,
        qc_h_qubits,
        k,
    )

    cs3.cx(4, 0)
    cs3.cx(4, 1)
    cs3.cx(4, 2)
    cs3.cx(4, 3)

    cs3.h(qc_h_qubits)
    cs3 = insert_error(
        cs3,
        qc_h_qubits,
        k,
    )

    # --------------------------------------------------------------
    # Fourth stabilizer-measurement subcircuit
    # --------------------------------------------------------------

    cs4 = QuantumCircuit(5)

    qc_h_qubits = [0, 3]

    cs4.h(qc_h_qubits)
    cs4 = insert_error(
        cs4,
        qc_h_qubits,
        k,
    )

    cs4.cx(4, 0)
    cs4.cx(4, 1)
    cs4.cx(4, 2)
    cs4.cx(4, 3)

    cs4.h(qc_h_qubits)
    cs4 = insert_error(
        cs4,
        qc_h_qubits,
        k,
    )

    # --------------------------------------------------------------
    # Compose the four stabilizer measurements
    # --------------------------------------------------------------

    qc = QuantumCircuit(
        9,
        name="$Error$ $Detection$",
    )

    qc.barrier()

    # Prepare the four syndrome ancillas in the |+> state.
    qc.h(range(5, 9))

    qc = insert_error(
        qc,
        range(5, 9),
        k,
    )

    qc.barrier()

    # Measure the first stabilizer using ancilla qubit 5.
    qc.compose(
        cs1,
        [4, 3, 2, 1, 5],
        inplace=True,
    )

    qc.barrier()

    # Measure the second stabilizer using ancilla qubit 6.
    qc.compose(
        cs2,
        [3, 2, 1, 0, 6],
        inplace=True,
    )

    qc.barrier()

    # Measure the third stabilizer using ancilla qubit 7.
    qc.compose(
        cs3,
        [0, 1, 2, 4, 7],
        inplace=True,
    )

    qc.barrier()

    # Measure the fourth stabilizer using ancilla qubit 8.
    qc.compose(
        cs4,
        [4, 3, 1, 0, 8],
        inplace=True,
    )

    qc.barrier()

    # Rotate the syndrome ancillas back to the computational basis.
    qc.h(range(5, 9))

    qc = insert_error(
        qc,
        range(5, 9),
        k,
    )

    qc.barrier()

    return qc


# ---------------------------------------------------------------------
# Syndrome decoding and recovery
# ---------------------------------------------------------------------

def decode_correction(qc, qbits, cbits):
    """
    Measure the syndrome and apply recovery for the [[5,1,3]] code.

    The four syndrome bits are decoded using a hard-decision lookup table.
    Each nonzero syndrome is mapped to one single-qubit Pauli correction.

    Syndrome qubits
    ----------------
    qbits[5], qbits[6], qbits[7], qbits[8]

    Syndrome classical bits
    -----------------------
    cbits[0], cbits[1], cbits[2], cbits[3]

    Syndrome ordering
    -----------------
    c0 c1 c2 c3

    Parameters
    ----------
    qc : QuantumCircuit
        Circuit containing the data and syndrome qubits.

    qbits : QuantumRegister
        Quantum register containing five data qubits followed by four
        syndrome ancillas.

    cbits : sequence of Clbit
        Four classical bits used to store the measured syndrome.

    Returns
    -------
    QuantumCircuit
        Circuit containing syndrome measurement, conditional recovery,
        and ancilla reset.
    """

    # Measure the four syndrome ancillas.
    for i in range(4):
        qc.measure(
            qbits[i + 5],
            cbits[i],
        )

    # -----------------------------------------------------------------
    # Syndrome lookup table
    #
    # 0000 corresponds to no detected error and therefore requires no
    # correction. Every other listed syndrome is mapped to a single-qubit
    # X, Y, or Z correction.
    # -----------------------------------------------------------------

    # First syndrome bit is 0: 0---
    with qc.if_test((cbits[0], 0)) as else_0:

        # First two syndrome bits are 00: 00--
        with qc.if_test((cbits[1], 0)) as else_1:

            # First three syndrome bits are 000: 000-
            with qc.if_test((cbits[2], 0)) as else_2:

                # Syndrome 0001: apply X to data qubit 4.
                with qc.if_test((cbits[3], 1)):
                    qc.x(qbits[4])

            # First three syndrome bits are 001: 001-
            with else_2:

                # Syndrome 0010: apply Z to data qubit 2.
                with qc.if_test((cbits[3], 0)) as else_3:
                    qc.z(qbits[2])

                # Syndrome 0011: apply X to data qubit 0.
                with else_3:
                    qc.x(qbits[0])

        # First two syndrome bits are 01: 01--
        with else_1:

            # First three syndrome bits are 010: 010-
            with qc.if_test((cbits[2], 0)) as else_4:

                # Syndrome 0100: apply Z to data qubit 0.
                with qc.if_test((cbits[3], 0)) as else_5:
                    qc.z(qbits[0])

                # Syndrome 0101: apply Z to data qubit 3.
                with else_5:
                    qc.z(qbits[3])

            # First three syndrome bits are 011: 011-
            with else_4:

                # Syndrome 0110: apply X to data qubit 1.
                with qc.if_test((cbits[3], 0)) as else_6:
                    qc.x(qbits[1])

                # Syndrome 0111: apply Y to data qubit 0.
                with else_6:
                    qc.y(qbits[0])

    # First syndrome bit is 1: 1---
    with else_0:

        # First two syndrome bits are 10: 10--
        with qc.if_test((cbits[1], 0)) as else_7:

            # First three syndrome bits are 100: 100-
            with qc.if_test((cbits[2], 0)) as else_8:

                # Syndrome 1000: apply X to data qubit 3.
                with qc.if_test((cbits[3], 0)) as else_9:
                    qc.x(qbits[3])

                # Syndrome 1001: apply Z to data qubit 1.
                with else_9:
                    qc.z(qbits[1])

            # First three syndrome bits are 101: 101-
            with else_8:

                # Syndrome 1010: apply Z to data qubit 4.
                with qc.if_test((cbits[3], 0)) as else_10:
                    qc.z(qbits[4])

                # Syndrome 1011: apply Y to data qubit 4.
                with else_10:
                    qc.y(qbits[4])

        # First two syndrome bits are 11: 11--
        with else_7:

            # First three syndrome bits are 110: 110-
            with qc.if_test((cbits[2], 0)) as else_11:

                # Syndrome 1100: apply X to data qubit 2.
                with qc.if_test((cbits[3], 0)) as else_12:
                    qc.x(qbits[2])

                # Syndrome 1101: apply Y to data qubit 3.
                with else_12:
                    qc.y(qbits[3])

            # First three syndrome bits are 111: 111-
            with else_11:

                # Syndrome 1110: apply Y to data qubit 2.
                with qc.if_test((cbits[3], 0)) as else_13:
                    qc.y(qbits[2])

                # Syndrome 1111: apply Y to data qubit 1.
                with else_13:
                    qc.y(qbits[1])

    # Reset the syndrome ancillas so that they may be reused.
    for i in range(4):
        qc.reset(qbits[i + 5])

    return qc

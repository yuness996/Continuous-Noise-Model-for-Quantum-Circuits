"""
Error detection and correction circuits for the [[5,1,3]] code.

Recommended versions:
    qiskit >= 2.0.0
    qiskit-aer compatible with qiskit >= 2.0.0

This file contains:
    - the ideal syndrome-extraction circuit;
    - the standard recovery table;
    - the transformed-basis recovery table.

The correction routines use Qiskit's dynamic-circuit control flow
through QuantumCircuit.if_test(...).
"""

from qiskit import QuantumCircuit


def encode():
    """Encoding circuit of [[5,1,3]]"""
    qc = QuantumCircuit(5,name='$Encoding$')
    qc.h([4,3,2,1])
    qc.z(0)
    qc.cx(4,0)
    qc.cx(3,0)
    qc.cx(2,0)
    qc.cx(1,0)
    qc.cz(4,0)
    qc.cz(4,3)
    qc.cz(2,1)
    qc.cz(3,2)
    qc.cz(1,0)
    return qc

def Logical_Hadamard():
    """Logic Hadamard circuit of [[5,1,3]]"""
    qc = QuantumCircuit(5, name = '  H_L  ')
    qc.h(range(5))
    qc.barrier()
    qc.swap(0, 1)
    qc.swap(1, 3)
    qc.swap(1, 4)
    return qc

# ---------------------------------------------------------------------
# Syndrome-extraction blocks
# ---------------------------------------------------------------------

def make_syndrome_block(gates):
    """
    Build a 5-qubit syndrome-extraction block.

    Parameters
    ----------
    gates : list[tuple[str, int, int]]
        Two-qubit gates written as

            ("cx", control, target)
            ("cz", control, target)

    Returns
    -------
    QuantumCircuit
        A 5-qubit block.
    """
    block = QuantumCircuit(5)

    for gate, control, target in gates:
        if gate == "cx":
            block.cx(control, target)

        elif gate == "cz":
            block.cz(control, target)

        else:
            raise ValueError(f"Unknown gate type: {gate}")

    return block


def syndrome_blocks_513():
    """
    Return the syndrome-extraction blocks for the [[5,1,3]] code.
    """
    block_1 = make_syndrome_block([
        ("cx", 4, 0),
        ("cz", 4, 1),
        ("cz", 4, 2),
        ("cx", 4, 3),
    ])

    block_3 = make_syndrome_block([
        ("cz", 4, 0),
        ("cz", 4, 1),
        ("cx", 4, 2),
        ("cx", 4, 3),
    ])

    block_4 = make_syndrome_block([
        ("cz", 4, 0),
        ("cx", 4, 1),
        ("cx", 4, 2),
        ("cz", 4, 3),
    ])

    return block_1, block_3, block_4
    
# ---------------------------------------------------------------------
# Error-detection circuit
# ---------------------------------------------------------------------

def detection(add_barriers=True):
    """
    Ideal error-detection circuit of the [[5,1,3]] code.

    Qubits 0--4 are data qubits.
    Qubits 5--8 are syndrome ancillas.

    Returns
    -------
    QuantumCircuit
        A 9-qubit syndrome-extraction circuit.
    """
    block_1, block_3, block_4 = syndrome_blocks_513()

    qc = QuantumCircuit(9, name="$Error$ $Detection$")

    if add_barriers:
        qc.barrier()

    qc.h(range(5, 9))

    if add_barriers:
        qc.barrier()

    qc.compose(block_1, [4, 3, 2, 1, 5], inplace=True)

    if add_barriers:
        qc.barrier()

    qc.compose(block_1, [3, 2, 1, 0, 6], inplace=True)

    if add_barriers:
        qc.barrier()

    qc.compose(block_3, [0, 1, 2, 4, 7], inplace=True)

    if add_barriers:
        qc.barrier()

    qc.compose(block_4, [4, 3, 1, 0, 8], inplace=True)

    if add_barriers:
        qc.barrier()

    qc.h(range(5, 9))

    if add_barriers:
        qc.barrier()

    return qc
    
# ---------------------------------------------------------------------
# Recovery tables
# ---------------------------------------------------------------------

def decode_correction(qc, qbits, cbits):
    """
    Dynamic post-processing error correction for the [[5,1,3]] code.

    Syndrome qubits:
        qbits[5], qbits[6], qbits[7], qbits[8]

    Syndrome classical bits:
        cbits[0], cbits[1], cbits[2], cbits[3]

    The syndrome order is c0 c1 c2 c3.
    """

    # Measure syndrome qubits.
    for i in range(4):
        qc.measure(qbits[i + 5], cbits[i])

    # Syndrome 0---
    with qc.if_test((cbits[0], 0)) as else_0:

        # Syndrome 00--
        with qc.if_test((cbits[1], 0)) as else_1:

            # Syndrome 000-
            with qc.if_test((cbits[2], 0)) as else_2:

                # 0001 -> X4
                with qc.if_test((cbits[3], 1)):
                    qc.x(qbits[4])

            # Syndrome 001-
            with else_2:

                # 0010 -> Z2
                with qc.if_test((cbits[3], 0)) as else_3:
                    qc.z(qbits[2])

                # 0011 -> X0
                with else_3:
                    qc.x(qbits[0])

        # Syndrome 01--
        with else_1:

            # Syndrome 010-
            with qc.if_test((cbits[2], 0)) as else_4:

                # 0100 -> Z0
                with qc.if_test((cbits[3], 0)) as else_5:
                    qc.z(qbits[0])

                # 0101 -> Z3
                with else_5:
                    qc.z(qbits[3])

            # Syndrome 011-
            with else_4:

                # 0110 -> X1
                with qc.if_test((cbits[3], 0)) as else_6:
                    qc.x(qbits[1])

                # 0111 -> Y0
                with else_6:
                    qc.y(qbits[0])

    # Syndrome 1---
    with else_0:

        # Syndrome 10--
        with qc.if_test((cbits[1], 0)) as else_7:

            # Syndrome 100-
            with qc.if_test((cbits[2], 0)) as else_8:

                # 1000 -> X3
                with qc.if_test((cbits[3], 0)) as else_9:
                    qc.x(qbits[3])

                # 1001 -> Z1
                with else_9:
                    qc.z(qbits[1])

            # Syndrome 101-
            with else_8:

                # 1010 -> Z4
                with qc.if_test((cbits[3], 0)) as else_10:
                    qc.z(qbits[4])

                # 1011 -> Y4
                with else_10:
                    qc.y(qbits[4])

        # Syndrome 11--
        with else_7:

            # Syndrome 110-
            with qc.if_test((cbits[2], 0)) as else_11:

                # 1100 -> X2
                with qc.if_test((cbits[3], 0)) as else_12:
                    qc.x(qbits[2])

                # 1101 -> Y3
                with else_12:
                    qc.y(qbits[3])

            # Syndrome 111-
            with else_11:

                # 1110 -> Y2
                with qc.if_test((cbits[3], 0)) as else_13:
                    qc.y(qbits[2])

                # 1111 -> Y1
                with else_13:
                    qc.y(qbits[1])

    # Reset syndrome qubits.
    for i in range(4):
        qc.reset(qbits[i + 5])

    return qc


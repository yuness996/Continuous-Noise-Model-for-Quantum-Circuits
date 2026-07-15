"""Logical post-processing and repeated simulation for the [[7,1,3]] code."""

from qiskit_aer import AerSimulator


# Computational-basis support of the logical |0_L> state.
L0 = [
    "0000000", "1010101", "0110011", "1100110",
    "0001111", "1011010", "0111100", "1101001",
]

# Computational-basis support of the logical |1_L> state.
L1 = [
    "1111111", "0101010", "1001100", "0011001",
    "1110000", "0100101", "1000011", "0010110",
]


def post_processing_decoding(counts):
    """
    Decode physical measurement outcomes into one- or two-qubit
    logical outcomes for the [[7,1,3]] code.
    """
    if not counts:
        raise ValueError("The count dictionary is empty.")

    key_length = len(next(iter(counts)))

    # One encoded logical qubit, with or without syndrome bits.
    if key_length < 14:
        logical_counts = {"0": 0, "1": 0}

        for bit_string, count in counts.items():
            data_state = bit_string[::-1][:7][::-1]

            if data_state in L0:
                logical_counts["0"] += count
            else:
                # Preserve the behavior of the original code: every state
                # outside L0 is classified as logical one.
                logical_counts["1"] += count

    # Two encoded logical qubits.
    elif key_length in (14, 20):
        logical_counts = {"00": 0, "01": 0, "10": 0, "11": 0}

        for bit_string, count in counts.items():
            if key_length == 20:
                # Extract the two seven-qubit data blocks while skipping
                # the syndrome bits present between them.
                state_qubit_1 = bit_string[::-1][:7][::-1]
                state_qubit_2 = bit_string[::-1][13:20][::-1]
            else:
                state_qubit_1 = bit_string[7:]
                state_qubit_2 = bit_string[:7]

            if state_qubit_1 in L0 and state_qubit_2 in L0:
                logical_counts["00"] += count
            elif state_qubit_1 in L1 and state_qubit_2 in L1:
                logical_counts["11"] += count
            elif (state_qubit_1 in L0) or (state_qubit_2 in L1):
                logical_counts["01"] += count
            elif (state_qubit_1 in L1) or (state_qubit_2 in L0):
                logical_counts["10"] += count
            else:
                # Preserve the original fallback rule for states that
                # cannot be assigned to either logical support.
                logical_counts["01"] += count // 2
                logical_counts["10"] += count - count // 2
    else:
        raise ValueError(f"Unsupported measurement-key length: {key_length}.")

    return logical_counts


def QuantumSimulator(qc, shots=100):
    """Execute a circuit on the ideal Aer simulator and return counts."""
    simulator = AerSimulator()
    result = simulator.run(qc, shots=int(shots)).result()
    return result.get_counts()


def Sum_dict(dict1, dict2):
    """Merge two count dictionaries by summing matching entries."""
    merged_dict = dict1.copy()

    for key, value in dict2.items():
        merged_dict[key] = merged_dict.get(key, 0) + value

    return merged_dict


def runcircuitnoisy(qc, shots):
    """
    Execute the circuit in repeated 100-shot batches and accumulate counts.

    The function name is retained for compatibility. No Aer noise model is
    added here; any noise must already be inserted directly into ``qc``.
    """
    counts = {}

    for _ in range(int(shots)):
        batch_counts = QuantumSimulator(qc, shots=100)
        counts = Sum_dict(counts, batch_counts)

    return counts


def runcircxtimes(qc, shots, times):
    """Repeat the full simulation and logically decode each result."""
    results = []

    for _ in range(int(times)):
        physical_counts = runcircuitnoisy(qc, int(shots))
        results.append(post_processing_decoding(physical_counts))

    return results
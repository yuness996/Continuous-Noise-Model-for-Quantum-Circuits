"""Grover search under a continuous coherent single-qubit noise model.

Random coherent error rotations are inserted after every one-qubit gate in the
initial state preparation, decomposed phase oracle, and decomposed diffuser.
The output metric is the average probability of not measuring the marked state.
"""

import random
from typing import Dict, Iterable, List, Optional

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import ZGate
from qiskit_aer import AerSimulator


# Each call to ``quantum_simulator`` uses this many measurement shots.
SHOTS_PER_CALL = 100


def r_error_gate(
    circuit: QuantumCircuit,
    theta: float,
    phi: float,
    qubit: int,
) -> QuantumCircuit:
    """Append the coherent error unitary to one qubit.

    The implemented gate sequence is

        Z R_x(2 theta) R_z(-2 phi) Z.

    Parameters
    ----------
    circuit:
        Circuit to which the error gates are appended.
    theta, phi:
        Random rotation angles defining the coherent error.
    qubit:
        Index of the affected qubit.
    """
    circuit.z(qubit)
    circuit.rx(2 * theta, qubit)
    circuit.rz(-2 * phi, qubit)
    circuit.z(qubit)
    return circuit


def insert_error(
    circuit: QuantumCircuit,
    positions: Iterable[int],
    noise_strength: float,
) -> QuantumCircuit:
    """Insert independent coherent errors on the selected qubits.

    Both angles are sampled from a zero-mean Gaussian distribution with
    standard deviation ``noise_strength / 3``.
    """
    for qubit in positions:
        theta = np.random.normal(loc=0.0, scale=noise_strength / 3)
        phi = np.random.normal(loc=0.0, scale=noise_strength / 3)
        r_error_gate(circuit, theta, phi, qubit)

    return circuit


def phase_oracle(n_qubits: int, noise_strength: float, target: str) -> QuantumCircuit:
    """Build a noisy phase oracle that marks ``target``.

    Qiskit uses little-endian qubit ordering, so the target bit string is
    reversed before mapping its bits to qubit indices. The ideal oracle is
    decomposed, then a coherent error is appended after each one-qubit gate.
    """
    if len(target) != n_qubits:
        raise ValueError("The target length must equal the number of qubits.")

    multi_controlled_z = ZGate().control(n_qubits - 1)
    ideal = QuantumCircuit(n_qubits)

    reversed_target = target[::-1]
    for qubit, bit in enumerate(reversed_target):
        if bit == "0":
            ideal.x(qubit)

    ideal.compose(multi_controlled_z, range(n_qubits), inplace=True)

    for qubit, bit in enumerate(reversed_target):
        if bit == "0":
            ideal.x(qubit)

    # Expose the elementary one-qubit gates on which noise is applied.
    ideal = ideal.decompose(reps=10)

    noisy = QuantumCircuit(n_qubits)
    for instruction in ideal.data:
        noisy.append(instruction.operation, instruction.qubits, instruction.clbits)

        if len(instruction.qubits) == 1:
            qubit = noisy.find_bit(instruction.qubits[0]).index
            insert_error(noisy, [qubit], noise_strength)

    return noisy


def diffuser(n_qubits: int, noise_strength: float) -> QuantumCircuit:
    """Build the noisy Grover diffusion operator."""
    multi_controlled_z = ZGate().control(n_qubits - 1)
    ideal = QuantumCircuit(n_qubits)

    ideal.h(range(n_qubits))
    ideal.x(range(n_qubits))
    ideal.compose(multi_controlled_z, range(n_qubits), inplace=True)
    ideal.x(range(n_qubits))
    ideal.h(range(n_qubits))

    # Decompose before deciding which physical gates receive noise.
    ideal = ideal.decompose(reps=10)

    noisy = QuantumCircuit(n_qubits)
    for instruction in ideal.data:
        noisy.append(instruction.operation, instruction.qubits, instruction.clbits)

        if len(instruction.qubits) == 1:
            qubit = noisy.find_bit(instruction.qubits[0]).index
            insert_error(noisy, [qubit], noise_strength)

    return noisy


def build_grover_circuit(
    n_qubits: int,
    target: str,
    noise_strength: float,
) -> QuantumCircuit:
    """Construct the complete noisy Grover-search circuit."""
    circuit = QuantumCircuit(n_qubits)

    # Prepare the uniform superposition and add preparation noise.
    circuit.h(range(n_qubits))
    insert_error(circuit, range(n_qubits), noise_strength)
    circuit.barrier()

    # Near-optimal integer number of Grover iterations for one marked state.
    repetitions = int(np.pi / 4 * np.sqrt(2**n_qubits) - 0.5)

    for _ in range(repetitions):
        circuit.compose(
            phase_oracle(n_qubits, noise_strength, target),
            range(n_qubits),
            inplace=True,
        )
        circuit.barrier()

        circuit.compose(
            diffuser(n_qubits, noise_strength),
            range(n_qubits),
            inplace=True,
        )
        circuit.barrier()

    circuit.measure_all()
    return circuit


def quantum_simulator(
    circuit: QuantumCircuit,
    shots: int = SHOTS_PER_CALL,
) -> Dict[str, int]:
    """Run a circuit on the ideal Aer shot-based simulator."""
    backend = AerSimulator()
    transpiled = transpile(circuit, backend=backend, optimization_level=3)
    result = backend.run(transpiled, shots=shots).result()
    return result.get_counts()


def merge_counts(first: Dict[str, int], second: Dict[str, int]) -> Dict[str, int]:
    """Return the element-wise sum of two count dictionaries."""
    merged = first.copy()
    for bitstring, count in second.items():
        merged[bitstring] = merged.get(bitstring, 0) + count
    return merged


def run_repeated_batches(
    circuit: QuantumCircuit,
    batches: int,
    shots_per_batch: int = SHOTS_PER_CALL,
) -> Dict[str, int]:
    """Run the same sampled noisy circuit several times and merge the counts.

    The coherent angles are fixed because they are sampled while ``circuit`` is
    constructed. This function repeats measurements of that same circuit; it
    does not generate a new noisy circuit for each batch.
    """
    counts: Dict[str, int] = {}

    for _ in range(batches):
        batch_counts = quantum_simulator(circuit, shots_per_batch)
        counts = merge_counts(counts, batch_counts)

    return counts


def sampled_targets(n_qubits: int) -> List[str]:
    """Generate ``n_qubits`` representative target strings.

    For each Hamming weight of zeros from 1 to ``n_qubits``, one target is
    sampled uniformly from the strings with that number of zeros. This is not
    the full set of ``2**n_qubits`` computational-basis targets.
    """
    targets = []

    for number_of_zeros in range(1, n_qubits + 1):
        bits = ["1"] * n_qubits
        zero_positions = random.sample(range(n_qubits), number_of_zeros)

        for position in zero_positions:
            bits[position] = "0"

        targets.append("".join(bits))

    return targets


def run_grover_noise_sweep(
    exponent_min: float,
    exponent_max: float,
    num_points: int,
    batches: int,
    n_qubits: int,
    seed: Optional[int] = None,
) -> List[float]:
    """Compute the average Grover failure probability over a noise sweep.

    Noise strengths are logarithmically spaced between
    ``10**exponent_min`` and ``10**exponent_max``. For each strength, the
    result is averaged over the representative targets returned by
    :func:`sampled_targets`.
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    noise_strengths = np.logspace(exponent_min, exponent_max, num=num_points)
    failure_probabilities = []
    targets = sampled_targets(n_qubits)

    for noise_strength in noise_strengths:
        failures = 0
        total_measurements = 0

        for target in targets:
            circuit = build_grover_circuit(n_qubits, target, noise_strength)
            counts = run_repeated_batches(circuit, batches)

            target_successes = counts.get(target, 0)
            target_shots = batches * SHOTS_PER_CALL

            failures += target_shots - target_successes
            total_measurements += target_shots

        failure_probabilities.append(failures / total_measurements)

    return failure_probabilities


def save_results(values: List[float], filename: str, total_shots: int) -> None:
    """Append one result row and its shot count to a text file."""
    with open(filename, "a", encoding="utf-8") as output_file:
        output_file.write(f"Shots = {total_shots}\n")
        output_file.write(" ".join(map(str, values)) + "\n")


if __name__ == "__main__":
    number_of_batches = 200
    results = run_grover_noise_sweep(
        exponent_min=-3,
        exponent_max=-1,
        num_points=20,
        batches=number_of_batches,
        n_qubits=6,
    )

    save_results(
        results,
        filename="GSA6qubits_cont.txt",
        total_shots=number_of_batches * SHOTS_PER_CALL,
    )
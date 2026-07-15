from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, transpile
import numpy as np
import sympy as sp

from PostProcessing import *
from Oracles import *


# ---------------------------------------------------------
# Analytical final sigmas
# ---------------------------------------------------------

def calculate_sigmas_analytical(arr, k):
    """
    Calculate final local Gaussian widths analytically.

    Output
    ------
    q_sigmas[i] = (sigma_theta_i, sigma_phi_i)

    Model
    -----
    - Each noisy one-qubit gate adds independent Gaussian noise.
    - Variances add:
          sigma_new^2 = sigma_old^2 + k^2
    - H swaps theta-like and phi-like widths.
    - X, Y, Z, S, Sdg add fresh local noise but do not change the local widths.
    - CX does not add fresh noise.
    - SWAP exchanges the widths of the two qubits.

    This is the local independent approximation.
    It does not keep CNOT-generated two-qubit error strings.
    """

    n_qubits = arr[0]
    depth = len(arr) - 1

    # Store variances, not sigmas.
    var_params = [(sp.Integer(0), sp.Integer(0)) for _ in range(n_qubits)]

    def add_gate_noise(pair):
        vt, vp = pair
        return vt + k**2, vp + k**2

    def hadamard_pair(pair):
        vt, vp = pair
        return vp, vt

    for i in range(1, depth + 1):
        gate_type, payload = arr[i]

        if gate_type == "h":
            qubit = payload

            # H swaps X-like and Z-like local components.
            var_params[qubit] = hadamard_pair(var_params[qubit])

            # Fresh noise after H.
            var_params[qubit] = add_gate_noise(var_params[qubit])

        elif gate_type in ("x", "y", "z", "s", "sdg"):
            qubit = payload

            # Local-width approximation: only add fresh noise.
            var_params[qubit] = add_gate_noise(var_params[qubit])

        elif gate_type == "cx":
            # No fresh noise after CX in this approximation.
            continue

        elif gate_type == "swap":
            q0, q1 = payload
            var_params[q0], var_params[q1] = var_params[q1], var_params[q0]

        else:
            raise ValueError(f"Unknown gate type '{gate_type}'")

    # Convert variances to sigmas.
    q_sigmas = [
        (sp.sqrt(vt), sp.sqrt(vp))
        for vt, vp in var_params
    ]

    return q_sigmas


# ---------------------------------------------------------
# Sampling and error gates
# ---------------------------------------------------------

def choose(kt, kp, loct=0.0, locp=0.0, wrap=True):
    """
    Sample one pair of final angles.

        theta ~ N(loct, kt^2)
        phi   ~ N(locp, kp^2)

    kt and kp are standard deviations.
    """

    theta = np.random.normal(loc=loct, scale=kt)
    phi = np.random.normal(loc=locp, scale=kp)

    if wrap:
        theta = (theta + np.pi) % (2 * np.pi) - np.pi
        phi = (phi + np.pi) % (2 * np.pi) - np.pi

    return theta, phi


def r_error_gate(qc, theta, phi, qubit):
    """
    Continuous rotation error:

        Z Rx(2 theta) Rz(-2 phi) Z
    """

    qc.z(qubit)
    qc.rx(2 * theta, qubit)
    qc.rz(-2 * phi, qubit)
    qc.z(qubit)

    return qc


def insert_error(qc, qubit, theta, phi):
    """
    Insert one sampled continuous error on one qubit.
    """

    return r_error_gate(qc, theta, phi, qubit)


def insert_error11(qc, pos, k):
    """
    Insert fresh gate-level errors on one or more qubits.

    This is kept for the explicit noisy detection circuit.
    """

    for p in pos:
        theta = np.random.normal(loc=0.0, scale=k)
        phi = np.random.normal(loc=0.0, scale=k)
        qc = r_error_gate(qc, theta, phi, p)

    return qc


def sample_angles_from_sigmas(q_sigmas, k_symbol, sig):
    """
    Sample final independent angles from analytical sigmas.

    q_sigmas[i] = (sigma_theta_i, sigma_phi_i)

    Returns
    -------
    sampled_angles[i] = (theta_i, phi_i)
    """

    sampled_angles = []

    for st_expr, sp_expr in q_sigmas:
        st = float(st_expr.subs(k_symbol, sig))
        sp_ = float(sp_expr.subs(k_symbol, sig))

        theta, phi = choose(st, sp_)
        sampled_angles.append((theta, phi))

    return sampled_angles


# ---------------------------------------------------------
# Convert Qiskit circuit to array
# ---------------------------------------------------------

def circuit_to_array(qc):
    """
    Convert a Qiskit circuit into the gate-list format used by the
    analytical approximation.

    Output
    ------
    ops : list
        Gate list in the form:
            [[gate_name, payload], ...]

    single_qubit_counts : list[int]
        single_qubit_counts[q] is the number of single-qubit gates
        applied to qubit q.

    Notes
    -----
    CZ is replaced by:
        H(target) - CX(control,target) - H(target)

    Therefore, each CZ adds two H gates to the target count.
    """

    flat = transpile(
        qc,
        basis_gates=["h", "x", "y", "z", "s", "sdg", "cx", "cz", "swap"],
        optimization_level=0,
    )

    ops = []
    single_qubit_counts = [0 for _ in range(flat.num_qubits)]

    for inst in flat.data:
        name = inst.operation.name.lower()
        qubits = [flat.qubits.index(q) for q in inst.qubits]

        # These do not affect the local-width propagation.
        if name in ("barrier", "measure", "delay"):
            continue

        if name == "cz":
            control, target = qubits

            ops.append(["h", target])
            single_qubit_counts[target] += 1

            ops.append(["cx", [control, target]])

            ops.append(["h", target])
            single_qubit_counts[target] += 1

        elif len(qubits) == 1:
            q = qubits[0]

            ops.append([name, q])
            single_qubit_counts[q] += 1

        else:
            ops.append([name, qubits])

    return ops, single_qubit_counts


# ---------------------------------------------------------
# [[7,1,3]] encoded circuit
# ---------------------------------------------------------

def Circuit(n_h):
    """
    Build only the 7-data-qubit encoded logical circuit.

    This circuit is used only to compute the final analytical sigmas.
    The noisy detection circuit is not included here. It is added later
    explicitly in single_encoded_qubit().
    """

    n_data = 7

    qq = QuantumRegister(n_data)
    cc = ClassicalRegister(n_data)
    qc = QuantumCircuit(qq, cc)

    qc.compose(encode(), range(n_data), inplace=True)

    for _ in range(n_h):
        qc.append(Logical_Hadamard(), range(n_data))

    return qc


# ---------------------------------------------------------
# [[7,1,3]] noisy detection circuit
# ---------------------------------------------------------

def detection_ed(k):
    """
    Error detection circuit of the [[7,1,3]] code.

    This keeps the same detection structure as Tests31.py.
    Fresh coherent noise is inserted explicitly after the Hadamard gates.
    """

    qc = QuantumCircuit(13, name="Detection")

    qc.h(range(7, 13))
    qc = insert_error11(qc, range(7, 13), k)

    qc.barrier()
    qc.cx(7, 0)
    qc.cx(7, 2)
    qc.cx(7, 4)
    qc.cx(7, 6)

    qc.barrier()
    qc.cx(8, 0)
    qc.cx(8, 1)
    qc.cx(8, 4)
    qc.cx(8, 5)

    qc.barrier()
    qc.cx(9, 0)
    qc.cx(9, 1)
    qc.cx(9, 2)
    qc.cx(9, 3)

    qc.barrier()
    qc.h([0, 2, 4, 6])
    qc = insert_error11(qc, [0, 2, 4, 6], k)
    qc.cx(10, 0)
    qc.cx(10, 2)
    qc.cx(10, 4)
    qc.cx(10, 6)
    qc.h([0, 2, 4, 6])
    qc = insert_error11(qc, [0, 2, 4, 6], k)

    qc.barrier()
    qc.h([0, 1, 4, 5])
    qc = insert_error11(qc, [0, 1, 4, 5], k)
    qc.cx(11, 0)
    qc.cx(11, 1)
    qc.cx(11, 4)
    qc.cx(11, 5)
    qc.h([0, 1, 4, 5])
    qc = insert_error11(qc, [0, 1, 4, 5], k)

    qc.barrier()
    qc.h([0, 1, 2, 3])
    qc = insert_error11(qc, [0, 1, 2, 3], k)
    qc.cx(12, 0)
    qc.cx(12, 1)
    qc.cx(12, 2)
    qc.cx(12, 3)
    qc.h([0, 1, 2, 3])
    qc = insert_error11(qc, [0, 1, 2, 3], k)

    qc.barrier()
    qc.h(range(7, 13))
    qc = insert_error11(qc, range(7, 13), k)

    return qc


# ---------------------------------------------------------
# Final-error approximation circuit
# ---------------------------------------------------------

def single_encoded_qubit(shots, times, n_h, correction, sampled_angles, sig):
    """
    Run the final-error approximation circuit for the [[7,1,3]] code.

    Workflow
    --------
    1. Prepare the encoded state.
    2. Apply the remaining ideal logical operation:
           H_L^n_h = I   if n_h is even,
           H_L     if n_h is odd.
    3. Insert one final effective coherent error on each data qubit.
    4. If correction=True, add the noisy detection circuit, then decode
       and correct.
    5. Measure the 7 data qubits.

    sampled_angles must contain one (theta, phi) pair for each data qubit.
    """

    n_data = 7
    n_ancilla = 6 if correction else 0
    n_qubits = n_data + n_ancilla

    if len(sampled_angles) != n_data:
        raise ValueError(
            f"sampled_angles must contain {n_data} data-qubit pairs, "
            f"but got {len(sampled_angles)}"
        )

    qq = QuantumRegister(n_qubits)
    cc = ClassicalRegister(n_qubits)
    qc = QuantumCircuit(qq, cc)

    # Encode one logical qubit into 7 data qubits.
    qc.compose(encode(), range(n_data), inplace=True)

    # Keep only the ideal logical action after moving the noise to the end.
    if n_h % 2 == 1:
        qc.append(h_L(), range(n_data))

    # Insert one final effective error on each data qubit.
    for q in range(n_data):
        theta, phi = sampled_angles[q]
        qc = insert_error(qc, qq[q], theta, phi)

    # Add noisy syndrome extraction and correction.
    if correction:
        qc.compose(detection_ed(k=sig), range(n_qubits), inplace=True)

        qc = decode_correction(
            qc,
            qq,
            [cc[7], cc[8], cc[9], cc[10], cc[11], cc[12]],
        )

    # Measure only the data qubits.
    for q in range(n_data):
        qc.measure(qq[q], cc[q])

    return runcircxtimes(qc, shots, times)


# ---------------------------------------------------------
# Main noise scan
# ---------------------------------------------------------

def Run1Qbt_AllNoise(a, b, num, shots, times, n_h, correction):
    """
    Run the [[7,1,3]] approximation model over a range of Gaussian widths.

    Steps
    -----
    1. Build the 7-data-qubit encoded logical circuit.
    2. Convert it to the gate-list format.
    3. Use calculate_sigmas_analytical to compute the final local widths.
    4. For each noise strength, sample one final error per data qubit.
    5. Run the final-error circuit.
    """

    n_data = 7

    # These are Gaussian standard deviations, not probabilities.
    sigma_values = (1 / 3) * np.logspace(a, b, num=num)

    k = sp.symbols("k", positive=True, real=True)

    # Compute final effective widths only for the encoded data block.
    ops, _ = circuit_to_array(Circuit(n_h))
    arr = [n_data] + ops

    q_sigmas = calculate_sigmas_analytical(arr, k)

    if len(q_sigmas) != n_data:
        raise ValueError(
            f"Expected {n_data} final sigma pairs, but got {len(q_sigmas)}"
        )

    results = []

    for sig in sigma_values:
        sampled_angles = sample_angles_from_sigmas(
            q_sigmas=q_sigmas,
            k_symbol=k,
            sig=sig,
        )

        result = single_encoded_qubit(
            shots=shots,
            times=times,
            n_h=n_h,
            correction=correction,
            sampled_angles=sampled_angles,
            sig=sig,
        )

        results.append(result)

    return results


# ---------------------------------------------------------
# Saving
# ---------------------------------------------------------

def save_table(table, filename, shots):
    """
    Save the raw result table in append mode.
    """

    with open(filename, "a") as file:
        file.write(f"Shots = {shots}\n")

        for row in table:
            if isinstance(row, (list, tuple, np.ndarray)):
                file.write(" ".join(map(str, row)) + "\n")
            else:
                file.write(str(row) + "\n")

        file.write("\n")


# ---------------------------------------------------------
# Run
# ---------------------------------------------------------

shots = 10
_times = 2
n_h = 2
correction = True

if correction:
    ec = "1"
else:
    ec = "0"

table = Run1Qbt_AllNoise(
    a=-6,
    b=-1,
    num=20,
    shots=shots,
    times=_times,
    n_h=n_h,
    correction=correction,
)

save_table(
    table,
    filename=f"Enc713Approx-m={str(n_h)}-ec={ec}_ed5.txt",
    shots=shots,
)

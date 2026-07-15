import numpy as np
from qiskit import QuantumCircuit
import random
def random_clifford_circuit_fixed_counts(
    n_qubits: int,
    n_cx: int,
    n_h: int,
    seed: int = None
):
    if n_h < n_qubits:
        raise ValueError("n_h must be >= n_qubits to ensure one single-qubit gate per qubit")

    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    circuit = [n_qubits]

    # 1. guarantee one H per qubit
    for q in range(n_qubits):
        circuit.append(("h", q))

    # 2. remaining H gates
    for _ in range(n_h - n_qubits):
        q = random.randrange(n_qubits)
        circuit.append(("h", q))

    # 3. CX gates
    for _ in range(n_cx):
        control, target = random.sample(range(n_qubits), 2)
        circuit.append(("cx", [control, target]))

    # 4. randomize order (excluding n_qubits header)
    circuit[1:] = random.sample(circuit[1:], len(circuit) - 1)

    return circuit

import sympy as sp
#def calculate_parameters_symbolic(arr, k, v_par):
def calculate_sigmas_analytical(arr, k):
    """
    Calculate final local Gaussian widths analytically.

    Output:
        q_params[i] = (sigma_theta_i, sigma_phi_i)

    Model:
        - errors are inserted after each single-qubit gate h or s
        - each fresh error adds variance k^2 to theta and phi
        - H swaps theta-like and phi-like errors
        - S is treated as adding fresh local noise only
        - CX does not add noise and is ignored in this local approximation
    """
    n_qubits = arr[0]
    depth = len(arr) - 1

    # Store variances first, not sigmas
    var_params = [(sp.Integer(0), sp.Integer(0)) for _ in range(n_qubits)]

    def add_gate_noise(pair):
        vt, vp = pair
        return (vt + k**2, vp + k**2)

    def hadamard_pair(pair):
        vt, vp = pair
        return (vp, vt)

    for i in range(1, depth + 1):
        gate_type, payload = arr[i]

        if gate_type == "h":
            qubit = payload

            # propagate old errors through H
            var_params[qubit] = hadamard_pair(var_params[qubit])

            # add fresh error after H
            var_params[qubit] = add_gate_noise(var_params[qubit])

        elif gate_type == "s":
            qubit = payload

            # local approximation: do not track X -> Y
            # only add fresh error after S
            var_params[qubit] = add_gate_noise(var_params[qubit])

        elif gate_type == "cx":
            # no fresh physical noise after CX in your model
            # local approximation ignores two-qubit propagated terms
            continue

        else:
            raise ValueError(f"Unknown gate type '{gate_type}'")

    # Convert variances to sigmas
    q_params = [
        (sp.sqrt(vt), sp.sqrt(vp))
        for vt, vp in var_params
    ]

    return q_params

def choose(kt, kp, loct=0.0, locp=0.0, wrap=True):
    """
    Sample one pair of angles:
        theta ~ N(loct, kt^2)
        phi   ~ N(locp, kp^2)

    kt and kp are standard deviations, not variances.
    """
    theta = np.random.normal(loc=loct, scale=kt)
    phi   = np.random.normal(loc=locp, scale=kp)

    if wrap:
        theta = (theta + np.pi) % (2*np.pi) - np.pi
        phi   = (phi   + np.pi) % (2*np.pi) - np.pi

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
    Insert one sampled error on one qubit.
    """
    return r_error_gate(qc, theta, phi, qubit)


def insert_error1(qc, qubits, sigma):
    """
    Insert independent fresh errors on selected qubits.
    This is used for the full noisy model after each noisy gate.
    """
    for q in qubits:
        theta, phi = choose(sigma, sigma)
        qc = r_error_gate(qc, theta, phi, q)

    return qc


def sample_angles_from_sigmas(q_sigmas, k_symbol, sig):
    """
    q_sigmas[i] = (sigma_theta_i, sigma_phi_i), possibly symbolic in k.

    Returns:
        sampled_angles[i] = (theta_i, phi_i)
    """
    sampled_angles = []

    for st_expr, sp_expr in q_sigmas:
        st = float(st_expr.subs(k_symbol, sig))
        sp = float(sp_expr.subs(k_symbol, sig))

        theta, phi = choose(st, sp)
        sampled_angles.append((theta, phi))

    return sampled_angles


def construct_random_circuit(arr, sig, sampled_angles):
    """
    Build three circuits:

    qc_ideal:
        ideal circuit.

    qc_model:
        full noisy model. Fresh errors are inserted after each single-qubit gate.

    qc_apprxm:
        approximation model. No errors during the circuit.
        Final sampled errors are inserted at the end.

    Assumption:
        CNOT does not receive fresh physical noise in this model.
    """
    n_qubits = arr[0]
    depth = len(arr) - 1

    if len(sampled_angles) != n_qubits:
        raise ValueError(
            f"sampled_angles must contain {n_qubits} pairs, "
            f"but got {len(sampled_angles)}"
        )

    qc_ideal = QuantumCircuit(n_qubits)
    qc_model = QuantumCircuit(n_qubits)
    qc_apprxm = QuantumCircuit(n_qubits)

    for i in range(1, depth + 1):
        gate_type, payload = arr[i]

        if gate_type == "h":
            qubit = payload

            qc_ideal.h(qubit)

            qc_model.h(qubit)
            qc_model = insert_error1(qc_model, [qubit], sig)

            qc_apprxm.h(qubit)

        elif gate_type == "s":
            qubit = payload

            qc_ideal.s(qubit)

            qc_model.s(qubit)
            qc_model = insert_error1(qc_model, [qubit], sig)

            qc_apprxm.s(qubit)

        elif gate_type == "cx":
            control, target = payload

            qc_ideal.cx(control, target)
            qc_model.cx(control, target)
            qc_apprxm.cx(control, target)

        else:
            raise ValueError(f"Unknown gate type: {gate_type}")

    # final effective errors for the approximation circuit
    for q in range(n_qubits):
        theta_q, phi_q = sampled_angles[q]
        qc_apprxm = insert_error(qc_apprxm, q, theta_q, phi_q)

    return qc_ideal, qc_model, qc_apprxm

from qiskit_aer import StatevectorSimulator
from qiskit.quantum_info import Statevector
def StateVector(qc):
    """
    Return the final statevector of a quantum circuit.

    Parameters
    ----------
    qc : QuantumCircuit
        The circuit to simulate.

    Returns
    -------
    statevector : qiskit.quantum_info.Statevector or array-like
        The final statevector produced by the circuit.

    Notes
    -----
    This uses Qiskit's StatevectorSimulator. It is suitable here because
    your circuits contain only unitary gates and no intermediate measurements.
    """
    simulator = StatevectorSimulator()
    result = simulator.run(qc).result()
    statevector = result.get_statevector()
    return statevector


def fidelity(psi, phi) -> float:
    """
    Compute the pure-state fidelity between two statevectors.

    For two pure states |psi> and |phi>, the fidelity is

        F = |<psi|phi>|^2.

    Parameters
    ----------
    psi, phi : Statevector or array-like
        The two statevectors to compare.

    Returns
    -------
    float
        Fidelity between 0 and 1.
    """

    # If Qiskit Statevector objects are given, extract the NumPy arrays.
    if isinstance(psi, Statevector):
        psi = psi.data

    if isinstance(phi, Statevector):
        phi = phi.data

    # Convert inputs to complex NumPy arrays.
    psi = np.asarray(psi, dtype=np.complex128)
    phi = np.asarray(phi, dtype=np.complex128)

    # Check that both objects are valid statevectors.
    if psi.ndim != 1 or phi.ndim != 1:
        raise ValueError("Inputs must be 1D statevectors.")

    if psi.shape != phi.shape:
        raise ValueError("Shape mismatch between the two statevectors.")

    # Normalize both states to avoid numerical errors.
    norm_psi = np.linalg.norm(psi)
    norm_phi = np.linalg.norm(phi)

    if norm_psi == 0 or norm_phi == 0:
        raise ValueError("Statevector norm cannot be zero.")

    psi = psi / norm_psi
    phi = phi / norm_phi

    # np.vdot conjugates the first argument, so this is <psi|phi>.
    return float(np.abs(np.vdot(psi, phi))**2)


def mean_ratio(a, b):
    """
    Compute the mean of element-wise ratios:

        mean(a_i / b_i)

    Parameters
    ----------
    a, b : list or array-like
        Lists of equal length. In your case, these can be infidelities from
        the approximation and the full model.

    Returns
    -------
    float
        Average ratio.

    Notes
    -----
    This is useful for computing

        mean(infidelity_approx / infidelity_full).

    Be careful if some values of b are very small, because the ratio can
    become unstable.
    """

    if len(a) != len(b):
        raise ValueError("Lists must have the same length.")

    ratios = []

    for x, y in zip(a, b):
        if y == 0:
            raise ZeroDivisionError("A denominator in b is zero.")

        ratios.append(x / y)

    return float(np.mean(ratios))


def Run_Circuit(
    n_qubits: int,
    n_cx: int,
    n_h: int,
    sigma: float = 1e-4,
    n_noise_samples: int = 100,
):
    """
    Compare the full noisy model with the approximation model.

    Full model:
        - Build the random Clifford circuit.
        - Insert fresh rotational errors after each noisy single-qubit gate.

    Approximation model:
        - Calculate final qubit-dependent sigmas analytically.
        - Sample one final error per qubit from these sigmas.
        - Insert these errors only at the end of the circuit.

    Returns
    -------
    ratio : float
        infidelity_approx / infidelity_full

    infidelity_full : float
        1 - average fidelity of the full noisy model.

    infidelity_approx : float
        1 - average fidelity of the approximation model.
    """

    # Random Clifford circuit
    # If your generator uses n_sqg instead of n_h, keep n_sqg=n_h.
    arr = random_clifford_circuit_fixed_counts(
        n_qubits=n_qubits,
        n_cx=n_cx,
        n_h=n_h,
        seed=None
    )

    # Symbolic noise scale
    k = sp.symbols("k", positive=True, real=True)

    # Analytical final sigmas for each qubit
    # q_sigmas[i] = (sigma_theta_i, sigma_phi_i)
    q_sigmas = calculate_sigmas_analytical(arr, k)

    Fid_model = 0.0
    Fid_apprxm = 0.0

    for _ in range(n_noise_samples):

        # Sample final approximation angles from the analytical sigmas
        sampled_angles = sample_angles_from_sigmas(
            q_sigmas=q_sigmas,
            k_symbol=k,
            sig=sigma
        )

        # Build ideal, full noisy, and approximation circuits
        qc_ideal, qc_model, qc_apprxm = construct_random_circuit(
            arr=arr,
            sig=sigma,
            sampled_angles=sampled_angles
        )

        # Statevectors
        ideal = StateVector(qc_ideal)
        model = StateVector(qc_model)
        apprxm = StateVector(qc_apprxm)

        # Accumulate fidelities
        Fid_model += fidelity(ideal, model)
        Fid_apprxm += fidelity(ideal, apprxm)

    # Average fidelities
    Fid_model /= n_noise_samples
    Fid_apprxm /= n_noise_samples

    # Infidelities
    infidelity_full = 1.0 - Fid_model
    infidelity_approx = 1.0 - Fid_apprxm

    if infidelity_full == 0:
        raise ZeroDivisionError("Full-model infidelity is zero, ratio is undefined.")

    ratio = infidelity_approx / infidelity_full

    return ratio, infidelity_full, infidelity_approx


def Run_different_circuits(
    n_qubits: int,
    sigma: float = 1e-4,
    n_noise_samples: int = 100
):
    """
    Run the approximation test on a 10x10 grid.

    Rows:
        n_cx = 10, 20, ..., 100

    Columns:
        n_h = 10, 20, ..., 100

    Returns
    -------
    ratios : list[list[float]]
        Grid of mean ratios:
            mean(infidelity_approx / infidelity_full)
    """

    ratios = [[0.0 for _ in range(10)] for _ in range(10)]

    for i_index, n_cx in enumerate(range(10, 110, 10)):
        for j_index, n_h in enumerate(range(10, 110, 10)):

            ratio, mean_full, mean_approx = Run_Circuit(
                n_qubits=n_qubits,
                n_cx=n_cx,
                n_h=n_h,
                sigma=sigma,
                n_noise_samples=n_noise_samples
            )

            ratios[i_index][j_index] = ratio

    return ratios
res = Run_different_circuits(10)

def save_table(table, filename):
    with open(filename, 'a') as file:
        for row in table:
            file.write(" ".join(map(str, row)) + "\n")

save_table(res, "RandomCircuits10qbts_new.txt")


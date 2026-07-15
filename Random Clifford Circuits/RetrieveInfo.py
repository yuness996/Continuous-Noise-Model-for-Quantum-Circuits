import random

import numpy as np
import sympy as sp

from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector
from qiskit_aer import StatevectorSimulator


# =====================================================================
# Random Clifford-circuit generation
# =====================================================================

def random_clifford_circuit_fixed_counts(
    n_qubits: int,
    n_cx: int,
    n_h: int,
    seed: int = None,
):
    """
    Generate a random Clifford circuit with fixed numbers of H and CX gates.

    At least one Hadamard gate is placed on each qubit. The remaining
    Hadamard and CX gates are generated randomly, then all gates are
    shuffled.

    The circuit is represented as a list rather than a QuantumCircuit:

        [
            n_qubits,
            ("h", qubit),
            ("cx", [control, target]),
            ...
        ]

    Parameters
    ----------
    n_qubits : int
        Number of qubits.
    n_cx : int
        Number of CNOT gates.
    n_h : int
        Total number of Hadamard gates.
    seed : int or None
        Optional random seed.

    Returns
    -------
    list
        Gate-list representation of the random circuit.
    """
    if n_qubits < 1:
        raise ValueError("n_qubits must be at least 1.")

    if n_cx < 0:
        raise ValueError("n_cx must be non-negative.")

    if n_h < n_qubits:
        raise ValueError(
            "n_h must be >= n_qubits to ensure at least one "
            "single-qubit gate per qubit."
        )

    if n_cx > 0 and n_qubits < 2:
        raise ValueError("At least two qubits are required for CX gates.")

    # Seed both random generators for reproducibility.
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    # The first entry stores the number of qubits.
    circuit = [n_qubits]

    # Ensure that every qubit receives at least one Hadamard gate.
    for qubit in range(n_qubits):
        circuit.append(("h", qubit))

    # Add the remaining Hadamard gates on random qubits.
    for _ in range(n_h - n_qubits):
        qubit = random.randrange(n_qubits)
        circuit.append(("h", qubit))

    # Add random CX gates with distinct control and target qubits.
    for _ in range(n_cx):
        control, target = random.sample(range(n_qubits), 2)
        circuit.append(("cx", [control, target]))

    # Randomize the gate order while keeping the qubit-count header first.
    circuit[1:] = random.sample(circuit[1:], len(circuit) - 1)

    return circuit


# =====================================================================
# Analytical propagation of local error widths
# =====================================================================

def calculate_sigmas_analytical(arr, k):
    """
    Calculate the final local Gaussian widths for each qubit.

    The calculation propagates local error variances through the circuit.
    Variances are tracked instead of standard deviations because
    independent Gaussian variances add directly.

    Model assumptions
    -----------------
    - A fresh error is inserted after every H or S gate.
    - Each fresh error adds k^2 to both local variances.
    - H exchanges the theta-like and phi-like components.
    - S only adds fresh noise in this local approximation.
    - CX adds no fresh noise.
    - Correlated error terms produced by CX propagation are ignored.

    Parameters
    ----------
    arr : list
        Gate-list representation of the circuit.
    k : sympy.Symbol
        Symbolic single-gate noise scale.

    Returns
    -------
    list[tuple]
        One pair ``(sigma_theta, sigma_phi)`` for each qubit.
    """
    n_qubits = arr[0]

    # Store variances:
    # var_params[q] = (variance_theta, variance_phi).
    var_params = [
        (sp.Integer(0), sp.Integer(0))
        for _ in range(n_qubits)
    ]

    def add_gate_noise(pair):
        """Add the variance of one fresh symmetric error."""
        var_theta, var_phi = pair

        return (
            var_theta + k**2,
            var_phi + k**2,
        )

    def propagate_through_hadamard(pair):
        """
        Propagate the two local error components through H.

        In this approximation, the theta-like and phi-like components
        are exchanged.
        """
        var_theta, var_phi = pair
        return var_phi, var_theta

    # Skip arr[0], which contains the number of qubits.
    for gate_type, payload in arr[1:]:

        if gate_type == "h":
            qubit = payload

            # Propagate errors already present before the H gate.
            var_params[qubit] = propagate_through_hadamard(
                var_params[qubit]
            )

            # Add the fresh error inserted after the H gate.
            var_params[qubit] = add_gate_noise(
                var_params[qubit]
            )

        elif gate_type == "s":
            qubit = payload

            # The local approximation does not explicitly track the
            # transformation X -> Y under S. Only fresh noise is added.
            var_params[qubit] = add_gate_noise(
                var_params[qubit]
            )

        elif gate_type == "cx":
            # CNOT gates receive no fresh physical noise.
            #
            # The full model can propagate existing errors through the
            # CNOT, but correlated multi-qubit terms are ignored here.
            continue

        else:
            raise ValueError(f"Unknown gate type '{gate_type}'.")

    # Convert the final variances into Gaussian standard deviations.
    q_sigmas = [
        (sp.sqrt(var_theta), sp.sqrt(var_phi))
        for var_theta, var_phi in var_params
    ]

    return q_sigmas


# =====================================================================
# Continuous-error sampling and insertion
# =====================================================================

def choose(
    sigma_theta,
    sigma_phi,
    mean_theta=0.0,
    mean_phi=0.0,
    wrap=True,
):
    """
    Sample one pair of continuous-error angles.

    The angles are sampled as

        theta ~ N(mean_theta, sigma_theta^2),
        phi   ~ N(mean_phi,   sigma_phi^2).

    Parameters
    ----------
    sigma_theta, sigma_phi : float
        Standard deviations of the two angles.
    mean_theta, mean_phi : float
        Gaussian means.
    wrap : bool
        If True, map the sampled angles to [-pi, pi).

    Returns
    -------
    tuple[float, float]
        Sampled ``(theta, phi)``.
    """
    theta = np.random.normal(
        loc=mean_theta,
        scale=sigma_theta,
    )

    phi = np.random.normal(
        loc=mean_phi,
        scale=sigma_phi,
    )

    if wrap:
        # Map large sampled angles to their equivalent values in [-pi, pi).
        theta = (theta + np.pi) % (2 * np.pi) - np.pi
        phi = (phi + np.pi) % (2 * np.pi) - np.pi

    return theta, phi


def r_error_gate(qc, theta, phi, qubit):
    """
    Apply the continuous rotational error to one qubit.

    The error unitary is decomposed as

        Z Rx(2 theta) Rz(-2 phi) Z.

    Parameters
    ----------
    qc : QuantumCircuit
        Circuit receiving the error.
    theta, phi : float
        Sampled error angles.
    qubit : int
        Qubit index.

    Returns
    -------
    QuantumCircuit
        Modified circuit.
    """
    qc.z(qubit)
    qc.rx(2 * theta, qubit)
    qc.rz(-2 * phi, qubit)
    qc.z(qubit)

    return qc


def insert_error(qc, qubit, theta, phi):
    """
    Insert one previously sampled error on one qubit.

    This helper is used by the approximation model, where one effective
    error is inserted at the end of the circuit.
    """
    return r_error_gate(
        qc=qc,
        theta=theta,
        phi=phi,
        qubit=qubit,
    )


def insert_error1(qc, qubits, sigma):
    """
    Insert independent fresh errors on selected qubits.

    This function is used by the full noisy model after every noisy
    single-qubit gate.

    Parameters
    ----------
    qc : QuantumCircuit
        Circuit receiving the errors.
    qubits : iterable[int]
        Qubits receiving independent errors.
    sigma : float
        Standard deviation of both sampled angles.

    Returns
    -------
    QuantumCircuit
        Modified circuit.
    """
    for qubit in qubits:
        theta, phi = choose(
            sigma_theta=sigma,
            sigma_phi=sigma,
        )

        qc = r_error_gate(
            qc=qc,
            theta=theta,
            phi=phi,
            qubit=qubit,
        )

    return qc


def sample_angles_from_sigmas(q_sigmas, k_symbol, sigma):
    """
    Sample final effective errors from symbolic local widths.

    Parameters
    ----------
    q_sigmas : list[tuple]
        Symbolic ``(sigma_theta, sigma_phi)`` for each qubit.
    k_symbol : sympy.Symbol
        Symbol appearing in the analytical expressions.
    sigma : float
        Numerical value substituted for the symbolic noise scale.

    Returns
    -------
    list[tuple]
        Sampled ``(theta, phi)`` for each qubit.
    """
    sampled_angles = []

    for sigma_theta_expr, sigma_phi_expr in q_sigmas:
        # Replace the symbolic gate-noise scale by its numerical value.
        sigma_theta = float(
            sigma_theta_expr.subs(k_symbol, sigma)
        )

        sigma_phi = float(
            sigma_phi_expr.subs(k_symbol, sigma)
        )

        theta, phi = choose(
            sigma_theta=sigma_theta,
            sigma_phi=sigma_phi,
        )

        sampled_angles.append((theta, phi))

    return sampled_angles


# =====================================================================
# Construction of the ideal, full, and approximate circuits
# =====================================================================

def construct_random_circuit(arr, sigma, sampled_angles):
    """
    Build the three circuits used in the comparison.

    Circuits
    --------
    qc_ideal
        Ideal circuit without errors.

    qc_model
        Full noisy model. A fresh continuous error is inserted after
        each single-qubit gate.

    qc_approx
        Approximate model. The ideal gates are applied first, then one
        effective sampled error is inserted on each qubit at the end.

    Assumptions
    -----------
    - H and S gates receive fresh continuous errors.
    - CX gates receive no fresh physical noise.
    - Errors can propagate naturally through later gates in the full model.
    - The approximation retains only independent final local errors.

    Parameters
    ----------
    arr : list
        Gate-list representation of the random circuit.
    sigma : float
        Standard deviation of each fresh gate-level error.
    sampled_angles : list[tuple]
        Final effective error angles for the approximate circuit.

    Returns
    -------
    tuple[QuantumCircuit, QuantumCircuit, QuantumCircuit]
        Ideal, full-model, and approximate-model circuits.
    """
    n_qubits = arr[0]

    if len(sampled_angles) != n_qubits:
        raise ValueError(
            f"sampled_angles must contain {n_qubits} pairs, "
            f"but got {len(sampled_angles)}."
        )

    qc_ideal = QuantumCircuit(n_qubits)
    qc_model = QuantumCircuit(n_qubits)
    qc_approx = QuantumCircuit(n_qubits)

    # Build the three circuits gate by gate.
    for gate_type, payload in arr[1:]:

        if gate_type == "h":
            qubit = payload

            # Ideal circuit.
            qc_ideal.h(qubit)

            # Full model: ideal H followed by a fresh error.
            qc_model.h(qubit)
            qc_model = insert_error1(
                qc=qc_model,
                qubits=[qubit],
                sigma=sigma,
            )

            # Approximate model: no error is inserted yet.
            qc_approx.h(qubit)

        elif gate_type == "s":
            qubit = payload

            # Ideal circuit.
            qc_ideal.s(qubit)

            # Full model: ideal S followed by a fresh error.
            qc_model.s(qubit)
            qc_model = insert_error1(
                qc=qc_model,
                qubits=[qubit],
                sigma=sigma,
            )

            # Approximate model: defer the error until the end.
            qc_approx.s(qubit)

        elif gate_type == "cx":
            control, target = payload

            # CX is ideal in all three circuits.
            qc_ideal.cx(control, target)
            qc_model.cx(control, target)
            qc_approx.cx(control, target)

        else:
            raise ValueError(f"Unknown gate type '{gate_type}'.")

    # Insert one final effective error on each qubit in the approximation.
    for qubit, (theta, phi) in enumerate(sampled_angles):
        qc_approx = insert_error(
            qc=qc_approx,
            qubit=qubit,
            theta=theta,
            phi=phi,
        )

    return qc_ideal, qc_model, qc_approx


# =====================================================================
# Statevector simulation and fidelity
# =====================================================================

def statevector(qc):
    """
    Return the final statevector of a unitary quantum circuit.

    The circuits used here contain no measurements, so statevector
    simulation can be used directly.
    """
    simulator = StatevectorSimulator()
    result = simulator.run(qc).result()

    return result.get_statevector()


def fidelity(psi, phi) -> float:
    """
    Compute the pure-state fidelity

        F = |<psi|phi>|^2.

    The inputs are normalized before the overlap is evaluated.

    Parameters
    ----------
    psi, phi : Statevector or array-like
        Statevectors to compare.

    Returns
    -------
    float
        Fidelity in the interval [0, 1], up to numerical precision.
    """
    # Extract arrays from Qiskit Statevector objects when needed.
    if isinstance(psi, Statevector):
        psi = psi.data

    if isinstance(phi, Statevector):
        phi = phi.data

    psi = np.asarray(psi, dtype=np.complex128)
    phi = np.asarray(phi, dtype=np.complex128)

    if psi.ndim != 1 or phi.ndim != 1:
        raise ValueError("Inputs must be one-dimensional statevectors.")

    if psi.shape != phi.shape:
        raise ValueError(
            "The two statevectors must have the same shape."
        )

    norm_psi = np.linalg.norm(psi)
    norm_phi = np.linalg.norm(phi)

    if norm_psi == 0 or norm_phi == 0:
        raise ValueError("A statevector norm cannot be zero.")

    # Normalize to reduce the effect of small numerical deviations.
    psi = psi / norm_psi
    phi = phi / norm_phi

    # np.vdot conjugates the first argument and evaluates <psi|phi>.
    overlap = np.vdot(psi, phi)

    return float(np.abs(overlap) ** 2)


def mean_ratio(a, b):
    """
    Compute the mean of the element-wise ratios ``a_i / b_i``.

    This helper is not used by the main simulation below, but can be used
    when several approximation and full-model infidelities are available.
    """
    if len(a) != len(b):
        raise ValueError("The two arrays must have the same length.")

    ratios = []

    for numerator, denominator in zip(a, b):
        if denominator == 0:
            raise ZeroDivisionError(
                "A denominator in the second array is zero."
            )

        ratios.append(numerator / denominator)

    return float(np.mean(ratios))


# =====================================================================
# Comparison for one random circuit
# =====================================================================

def run_circuit(
    n_qubits: int,
    n_cx: int,
    n_h: int,
    sigma: float = 1e-4,
    n_noise_samples: int = 100,
):
    """
    Compare the full and approximate error models for one random circuit.

    One random Clifford circuit is generated and kept fixed. Several
    independent noise samples are then applied to this same ideal circuit.

    Full model
    ----------
    A fresh rotational error is sampled after every noisy single-qubit gate.

    Approximate model
    -----------------
    The final local widths are calculated analytically. One effective
    rotational error is then sampled and applied at the end of each qubit.

    Parameters
    ----------
    n_qubits : int
        Number of qubits.
    n_cx : int
        Number of CNOT gates.
    n_h : int
        Number of Hadamard gates.
    sigma : float
        Standard deviation of each fresh gate-level error.
    n_noise_samples : int
        Number of independent noise realizations for this circuit.

    Returns
    -------
    ratio : float
        ``infidelity_approx / infidelity_full``.
    infidelity_full : float
        Mean infidelity of the full noisy model.
    infidelity_approx : float
        Mean infidelity of the approximate model.
    """
    if sigma < 0:
        raise ValueError("sigma must be non-negative.")

    if n_noise_samples < 1:
        raise ValueError("n_noise_samples must be at least 1.")

    # Generate one random ideal circuit structure.
    circuit_description = random_clifford_circuit_fixed_counts(
        n_qubits=n_qubits,
        n_cx=n_cx,
        n_h=n_h,
        seed=None,
    )

    # Symbolic gate-level noise scale used in the analytical calculation.
    k = sp.symbols("k", positive=True, real=True)

    # Final local effective widths:
    # q_sigmas[q] = (sigma_theta_q, sigma_phi_q).
    q_sigmas = calculate_sigmas_analytical(
        circuit_description,
        k,
    )

    fidelity_full_sum = 0.0
    fidelity_approx_sum = 0.0

    for _ in range(n_noise_samples):
        # Sample one effective final error for every qubit.
        sampled_angles = sample_angles_from_sigmas(
            q_sigmas=q_sigmas,
            k_symbol=k,
            sigma=sigma,
        )

        # Build the same ideal gate sequence under the two noise models.
        qc_ideal, qc_full, qc_approx = construct_random_circuit(
            arr=circuit_description,
            sigma=sigma,
            sampled_angles=sampled_angles,
        )

        # Simulate the three final pure states.
        psi_ideal = statevector(qc_ideal)
        psi_full = statevector(qc_full)
        psi_approx = statevector(qc_approx)

        # Compare both noisy states with the same ideal final state.
        fidelity_full_sum += fidelity(
            psi_ideal,
            psi_full,
        )

        fidelity_approx_sum += fidelity(
            psi_ideal,
            psi_approx,
        )

    # Average over independent noise realizations.
    mean_fidelity_full = fidelity_full_sum / n_noise_samples
    mean_fidelity_approx = fidelity_approx_sum / n_noise_samples

    # Convert fidelities into infidelities.
    infidelity_full = 1.0 - mean_fidelity_full
    infidelity_approx = 1.0 - mean_fidelity_approx

    if np.isclose(infidelity_full, 0.0):
        raise ZeroDivisionError(
            "The full-model infidelity is zero, so the ratio is undefined."
        )

    ratio = infidelity_approx / infidelity_full

    return ratio, infidelity_full, infidelity_approx


# =====================================================================
# Sweep over different circuit sizes
# =====================================================================

def run_different_circuits(
    n_qubits: int,
    sigma: float = 1e-4,
    n_noise_samples: int = 100,
):
    """
    Compare the models on a 10 x 10 grid of random circuit parameters.

    Rows correspond to:

        n_cx = 10, 20, ..., 100.

    Columns correspond to:

        n_h = 10, 20, ..., 100.

    A new random circuit is generated for every ``(n_cx, n_h)`` pair.

    Parameters
    ----------
    n_qubits : int
        Number of qubits in every circuit.
    sigma : float
        Gate-level continuous-noise standard deviation.
    n_noise_samples : int
        Number of noise realizations per random circuit.

    Returns
    -------
    list[list[float]]
        Grid containing

            infidelity_approx / infidelity_full

        for each circuit configuration.
    """
    ratios = [
        [0.0 for _ in range(10)]
        for _ in range(10)
    ]

    # Vary the number of CNOT gates along the rows.
    for row_index, n_cx in enumerate(range(10, 110, 10)):

        # Vary the number of Hadamard gates along the columns.
        for column_index, n_h in enumerate(range(10, 110, 10)):
            ratio, infidelity_full, infidelity_approx = run_circuit(
                n_qubits=n_qubits,
                n_cx=n_cx,
                n_h=n_h,
                sigma=sigma,
                n_noise_samples=n_noise_samples,
            )

            ratios[row_index][column_index] = ratio

    return ratios


# =====================================================================
# Saving results
# =====================================================================

def save_table(table, filename):
    """
    Append a two-dimensional numerical table to a text file.

    Append mode is used so that results from repeated or parallel jobs can
    be stored without deleting earlier data.

    Each row is written as space-separated numerical values.
    """
    with open(filename, "a", encoding="utf-8") as file:
        for row in table:
            file.write(" ".join(map(str, row)) + "\n")


# =====================================================================
# Main execution
# =====================================================================

if __name__ == "__main__":
    # Evaluate a 10 x 10 grid for ten-qubit random circuits.
    #
    # With the default arguments, each grid point is averaged over
    # 100 independent noise realizations.
    results = run_different_circuits(
        n_qubits=10,
        sigma=1e-4,
        n_noise_samples=100,
    )

    # Append the resulting matrix to the output file.
    save_table(
        results,
        "RandomCircuits10qbts_new.txt",
    )

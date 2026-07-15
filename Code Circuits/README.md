# Continuous Coherent Noise for the [[7,1,3]] Steane Code

This repository contains the Qiskit code used to build and test noisy quantum-error-correction circuits based on the `[[7,1,3]]` Steane code. Continuous coherent errors are modeled as small random single-qubit rotations and inserted after selected single-qubit gates.

The code includes the Steane-code encoding circuit, logical gates, syndrome extraction, hard-decision recovery, logical post-processing, and scripts for generating circuits over a range of noise strengths.

## Repository structure

```text
.
├── Oracles.py          # Circuit blocks(encoding, logical circuits, stabilizer measure), and recovery(decoding/correction)
├── PostProcessing.py   # Circuit execution and logical state retrieval
└── Tests.py            # Construction of circuits over a noise grid
```

## Files

### `Oracles.py`

Defines the quantum circuits and the continuous error model:

- `r_error_gate(qc, theta, phi, qubit)` applies the sampled coherent error rotation
  
  ```text
  Z Rx(2θ) Rz(-2φ) Z
  ```

- `insert_error(qc, pos, k)` samples independent Gaussian angles and inserts an error on each selected qubit.
- `encode(k)` builds the noisy encoding circuit for one `[[7,1,3]]` code block.
- `detection(k)` builds the six-ancilla syndrome-extraction circuit.
- `h_L(k)` applies the transversal logical Hadamard.
- `cnot()` applies a transversal logical CNOT between two encoded blocks.
- `decode_correction(qc, qbits, cbits)` measures the six syndrome ancillas and applies the corresponding single-qubit Pauli correction.

The sampled angles follow

```text
θ, φ ~ N(0, (k/3)²),
```

where `k` controls the noise strength.

### `PostProcessing.py`

Contains the simulation and logical decoding tools:

- `post_processing_decoding(counts)` maps physical measurement outcomes to logical outcomes.
- `QuantumSimulator(qc, shots)` runs a circuit on `AerSimulator`.
- `Sum_dict(dict1, dict2)` combines count dictionaries.
- `runcircuitnoisy(qc, shots)` repeats the simulation in 100-shot batches and accumulates the counts.
- `runcircxtimes(qc, shots, times)` repeats the full simulation and decodes each result.

Despite its name, `runcircuitnoisy` does not add an Aer noise model. The noise is already inserted directly into the circuit by the functions in `Oracles.py`.

### `Tests.py`

Builds the encoded circuit for a logarithmic range of coherent-noise strengths.

The current version:

1. encodes one logical qubit,
2. optionally applies one syndrome-extraction and correction cycle,
3. measures the seven data qubits,
4. returns the generated `QuantumCircuit` objects.

## Requirements

- Python 3.10 or later
- NumPy
- Qiskit
- Qiskit Aer

Install the required packages with

```bash
pip install numpy qiskit qiskit-aer
```

## Basic use

Build one encoded circuit without error correction:

```python
from Tests import single_encoded_qubit

qc = single_encoded_qubit(
    k=1e-3,
    shots=1000,
    times=10,
    correction=False,
)

print(qc)
```

Build one encoded circuit with syndrome extraction and correction:

```python
qc = single_encoded_qubit(
    k=1e-3,
    shots=1000,
    times=10,
    correction=True,
)
```

Run and decode the circuit:

```python
from PostProcessing import runcircxtimes

results = runcircxtimes(
    qc,
    shots=10,
    times=5,
)

print(results)
```

Here, `shots=10` means ten batches of 100 simulator shots, giving 1000 shots per repetition.

## Noise sweep

Generate circuits for 20 logarithmically spaced noise strengths between `10^-6` and `10^-1`:

```python
from Tests import Run1Qbt_AllNoise

circuits = Run1Qbt_AllNoise(
    a=-6,
    b=-1,
    num=20,
    shots=1000,
    times=10,
    correction=True,
)
```

The values of `shots` and `times` are currently retained for interface compatibility, but `Run1Qbt_AllNoise` returns circuits rather than simulation results.

## Important implementation notes

- Noise is added only after selected single-qubit gates.
- CNOT gates are treated as noiseless, but they propagate errors already present on the qubits.
- The logical Hadamard is transversal for the Steane code.
- The recovery step uses a fixed hard-decision syndrome lookup table.
- The post-processing function preserves the behavior of the original code: for one logical qubit, every measured state outside the support of `|0_L>` is classified as logical `1`. For stricter decoding, states outside both `L0` and `L1` should be tracked separately.
- In `h_L(k)`, coherent errors are currently inserted only on qubits `0, 1, 2, 3`, although Hadamard gates are applied to all seven qubits. This should be checked if noise is intended after every physical Hadamard.
- `Tests.py` currently builds circuits but does not execute them before saving. Circuit objects should be simulated first or saved using QPY/QASM rather than written as numerical table rows.

## Reproducibility

The noise angles are sampled randomly with NumPy. Set a seed before constructing the circuit to reproduce a given noisy instance:

```python
import numpy as np

np.random.seed(1234)
```

## Citation

This code accompanies the work:

> Y. El Kaderi, A. Honecker, and I. Andriyanova, *Continuous Noise Model for Quantum Circuits*.

Please cite the corresponding article when using this repository in academic work.

## License

Add a license file before public release. The MIT license is a common choice for research code intended for reuse.

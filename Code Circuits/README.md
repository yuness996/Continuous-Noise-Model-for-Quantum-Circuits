# Quantum-Circuit Simulation Codes

This directory contains the circuit-level simulations used to compare continuous coherent noise, its analytical approximation, and a Pauli depolarizing channel on the `[[5,1,3]]` and `[[7,1,3]]` quantum error-correcting codes.

Several files have the same names and similar structures. This is intentional. Each model and each code has its own self-contained set of:

* `Oracles.py`: circuit definitions;
* `PostProcessing.py`: simulation and logical decoding;
* `Tests.py`: parameter sweep and raw-data generation.

Do not mix files from different folders. Always run `Tests.py` from the folder containing the matching `Oracles.py` and `PostProcessing.py`.

---

## Which folder should I use?

| Goal                                                                        | Folder                                            |
| --------------------------------------------------------------------------- | ------------------------------------------------- |
| Simulate coherent errors explicitly after circuit gates                     | `Continuous Error Model/Full Model/`              |
| Use the faster local analytical approximation of coherent-error propagation | `Continuous Error Model/Approximate Error Model/` |
| Simulate a symmetric discrete Pauli channel                                 | `Pauli Depolarizing Channel/`                     |

Each model contains separate implementations for the `[[5,1,3]]` and `[[7,1,3]]` codes.

---

## Repository structure

```text
Code Circuits/
│
├── Continuous Error Model/
│   │
│   ├── Full Model/
│   │   ├── [[5,1,3]]ContinuousError/
│   │   │   ├── Oracles.py
│   │   │   ├── PostProcessing.py
│   │   │   └── Tests.py
│   │   │
│   │   └── [[7,1,3]]ContinuousError/
│   │       ├── Oracles.py
│   │       ├── PostProcessing.py
│   │       └── Tests.py
│   │
│   └── Approximate Error Model/
│       ├── [[5,1,3]] (continuous error) (model)/
│       │   ├── Oracles.py
│       │   ├── PostProcessing.py
│       │   └── Tests.py
│       │
│       └── [[7,1,3]] (continuous error) (model)/
│           ├── Oracles.py
│           ├── PostProcessing.py
│           └── Tests.py
│
├── Pauli Depolarizing Channel/
│   ├── [[5,1,3]] Pauli Errors/
│   │   ├── Oracles.py
│   │   ├── PostProcessing.py
│   │   └── Tests.py
│   │
│   └── [[7,1,3]] Pauli Errors/
│       ├── Oracles.py
│       ├── PostProcessing.py
│       └── Tests.py
│
├── RetrieveInformation.py
└── README.md
```

Folders named `__pycache__` may also appear. They are generated automatically by Python and are not part of the simulation workflow.

---

## Common role of the three Python files

### `Oracles.py`

This file defines the circuit blocks for one specific code and noise model. Depending on the folder, it contains:

* the encoding circuit;
* the logical Hadamard circuit;
* the syndrome-extraction circuit;
* the syndrome lookup table;
* conditional Pauli recovery;
* the continuous coherent-error gate, when the full continuous model is used.

The `[[5,1,3]]` implementation uses five data qubits and four syndrome ancillas when error correction is enabled.

The `[[7,1,3]]` Steane-code implementation uses seven data qubits and six syndrome ancillas when error correction is enabled.

### `PostProcessing.py`

This file runs the circuit and maps physical measurement results to logical outcomes.

It contains:

* the computational-basis supports of the logical states;
* conversion of physical counts into logical counts such as `{'0': ..., '1': ...}`;
* repeated circuit execution;
* accumulation of count dictionaries;
* the Qiskit Aer Pauli noise model in the Pauli folders.

In the continuous-model folders, the Aer simulator itself is ideal. The noise has already been inserted as gates in the circuit.

In the Pauli folders, `PostProcessing.py` creates an Aer `NoiseModel` and applies the Pauli channel during simulation.

### `Tests.py`

This is the main experiment file. It:

1. builds the encoded circuit;
2. applies the chosen number of logical Hadamard circuits;
3. optionally adds syndrome extraction and recovery;
4. scans a logarithmic range of noise strengths;
5. repeats the simulation;
6. saves the raw logical-count dictionaries to a text file.

This is normally the file that users should edit and execute.

---

## Noise models

### 1. Full continuous coherent-noise model

Location:

```text
Continuous Error Model/Full Model/
```

This model inserts coherent single-qubit rotations directly into the circuit. The error gate is

```text
Z Rx(2 theta) Rz(-2 phi) Z
```

where the angles are independently sampled from Gaussian distributions. In the current implementation,

```text
theta, phi ~ N(0, (k/3)^2)
```

and `k` controls the physical noise strength.

Fresh noise is inserted after selected single-qubit gates. Two-qubit gates do not receive a new error directly, although they can propagate errors already present on the qubits.

This is the most direct implementation of the continuous model. It explicitly samples the errors during circuit construction and is therefore more costly when many independent circuit instances are required.

---

### 2. Approximate continuous-error model

Location:

```text
Continuous Error Model/Approximate Error Model/
```

This version replaces the many gate-level coherent errors by final effective errors.

The main steps are:

1. build the ideal encoded circuit;
2. convert the circuit into a gate list;
3. track a pair of local variances for each qubit;
4. add the variance produced by noisy single-qubit gates;
5. propagate the local widths through the supported gates;
6. sample one final coherent rotation for each data qubit;
7. run the shorter final-error circuit.

The local propagation rules used in `Tests.py` include:

* a noisy single-qubit gate adds fresh variance;
* a Hadamard exchanges the two local error components;
* a SWAP exchanges the parameters of the two qubits;
* a CNOT does not add fresh noise;
* CNOT-generated multi-qubit correlations are not retained.

The approximation is faster because it avoids inserting and sampling every coherent error in the full circuit. It is a local independent approximation and should not be interpreted as an exact replacement for the full model.

When error correction is enabled, the effective accumulated error is applied to the data block, while the syndrome-extraction circuit contains its own explicitly inserted coherent errors.

---

### 3. Pauli depolarizing channel

Location:

```text
Pauli Depolarizing Channel/
```

These scripts use ideal circuit blocks and add a symmetric Pauli channel through Qiskit Aer.

For a total physical error probability `p`, the code assigns

```text
P(X) = p/3
P(Y) = p/3
P(Z) = p/3
P(I) = 1 - p
```

The current noise model is attached to the single-qubit basis gates `u1`, `u2`, and `u3`. It does not add fresh noise directly to CNOT gates.

This model is used for comparison with the continuous coherent-noise simulations.

---

## Installation

Python 3.10 or later is recommended.

Create a virtual environment and install the required packages:

```bash
python -m venv .venv
source .venv/bin/activate
pip install numpy sympy qiskit qiskit-aer
```

`SymPy` is required by the approximate continuous model.

The recovery circuits use Qiskit's dynamic control flow through `QuantumCircuit.if_test(...)`. Use compatible recent versions of Qiskit and Qiskit Aer.

---

## Running a simulation

First move to the exact experiment folder. This is important because the scripts use local imports such as:

```python
from Oracles import *
from PostProcessing import *
```

For example, to run the full continuous model for the `[[5,1,3]]` code:

```bash
cd "Code Circuits/Continuous Error Model/Full Model/[[5,1,3]]ContinuousError"
python Tests.py
```

To run the Pauli model for the `[[7,1,3]]` code:

```bash
cd "Code Circuits/Pauli Depolarizing Channel/[[7,1,3]] Pauli Errors"
python Tests.py
```

Do not run a `Tests.py` file while another experiment folder appears earlier in `PYTHONPATH`. Otherwise, Python may import the wrong `Oracles.py` or `PostProcessing.py`.

---

## Main simulation parameters

The parameters are set near the end of each `Tests.py` file.

A typical call has the form:

```python
results = Run_AllNoise(
    a=-6,
    b=-1,
    num=20,
    shots=1000,
    times=10,
    correction=True,
    n_h=10,
)
```

or, in the continuous-model folders:

```python
results = Run1Qbt_AllNoise(
    a=-6,
    b=-1,
    num=20,
    shots=1000,
    times=10,
    n_h=10,
    correction=True,
)
```

The parameters mean:

| Parameter    | Meaning                                                                   |
| ------------ | ------------------------------------------------------------------------- |
| `a`          | Lower base-10 exponent of the noise range                                 |
| `b`          | Upper base-10 exponent of the noise range                                 |
| `num`        | Number of logarithmically spaced noise values                             |
| `n_h`        | Number of logical Hadamard circuits                                       |
| `correction` | Enables syndrome extraction and recovery when `True`                      |
| `times`      | Number of repeated logical-count dictionaries stored for each noise value |
| `shots`      | Sampling parameter; its exact meaning depends on the model                |

For example,

```python
np.logspace(-6, -1, num=20)
```

creates 20 noise values between `10^-6` and `10^-1`.

---

## Important difference in the meaning of `shots`

The continuous and Pauli scripts do not currently use `shots` in the same way.

### Continuous full and approximate models

In the continuous `PostProcessing.py` files, `shots` is the number of batches, and each batch contains 100 Aer shots.

Therefore,

```text
effective shots per repeated dictionary = 100 x shots
```

For example:

```text
shots = 1000
```

produces `100000` physical measurement shots for each dictionary stored in the raw output.

### Pauli model

In the Pauli folders, `shots` is passed directly to Aer.

Therefore,

```text
effective shots per repeated dictionary = shots
```

Keep this difference in mind when comparing the raw count values or the total run cost.

---

## Why results are saved as raw text files

The simulations were run as independent parallel jobs because repeated noisy circuit sampling is expensive.

Each job computes a complete sweep over the noise grid and saves its raw logical-count dictionaries to a text file. Saving the raw data serves several purposes:

* each parallel job can finish independently;
* partial results are not lost if another job fails;
* several jobs can later be merged;
* error bars can be recomputed without rerunning the quantum circuits;
* the raw logical counts remain available for checks.

The simulation scripts use append mode:

```python
with open(filename, "a") as file:
```

This means that rerunning a script adds a new result block instead of replacing the old data.

Before starting a new experiment, either delete the old output file or choose a new filename. Otherwise, old and new jobs will be merged during post-processing.

---

## Recommended parallel workflow

It is safer for each parallel task to write to a different file. For example:

```text
Enc513Full-m=10-ec=1-job_001.txt
Enc513Full-m=10-ec=1-job_002.txt
Enc513Full-m=10-ec=1-job_003.txt
```

Here:

* `513` identifies the `[[5,1,3]]` code;
* `Full` identifies the full continuous model;
* `m=10` means ten logical Hadamard circuits;
* `ec=1` means error correction is enabled;
* `job_001` identifies the parallel task.

After all jobs finish, concatenate the files in a fixed order:

```bash
cat Enc513Full-m=10-ec=1-job_*.txt > Enc513Full-m=10-ec=1-combined.txt
```

Using separate files avoids the risk of several cluster tasks writing to the same file at the same time.

All files being merged must use the same:

* code;
* noise model;
* values of `a`, `b`, and `num`;
* number of logical gates;
* error-correction setting;
* number of repetitions;
* shot convention.

---

## Raw output format

A raw result file contains one row per noise value. Each row contains several logical-count dictionaries, normally one dictionary per repeated run.

A block may look like:

```text
Shots = 1000
{'0': 99990, '1': 10} {'0': 99987, '1': 13} ... {'0': 99991, '1': 9}
{'0': 99980, '1': 20} {'0': 99983, '1': 17} ... {'0': 99978, '1': 22}
...
```

The `Shots = ...` header is present in some scripts and absent in others. `RetrieveInformation.py` can skip these header lines.

With `num=20`, one complete job normally contributes 20 data rows. If several jobs are concatenated, the file contains several consecutive 20-row blocks.

These files contain raw data. They are not the final logical-error curves.

---

## Retrieving the final curve

The top-level script

```text
RetrieveInformation.py
```

reads the raw result blocks and returns the final logical-error probabilities and statistical error bars.

It performs the following steps:

1. ignores empty lines and `Shots = ...` headers;
2. extracts all dictionaries from each row;
3. groups every `sweep_size` rows into one parallel job;
4. merges matching noise points from all jobs by adding their count dictionaries;
5. computes the logical error probability for every repeated dictionary;
6. averages the probabilities;
7. computes the standard error of the mean.

For one encoded logical qubit, the logical error probability is

```text
P_L = N_1 / (N_0 + N_1)
```

where `N_0` and `N_1` are the decoded logical counts.

The header value is not used for normalization. The code uses the sum of the values inside each dictionary.

---

## Using `RetrieveInformation.py`

The simplest method is to edit the example at the bottom of the file:

```python
fn = "Enc513Full-m=10-ec=1-combined.txt"

curve, error_bar = load_error_curve(
    filename=fn,
    sweep_size=20,
    times=10,
    prepend_zero=True,
)

print("curve = [" + ",".join(map(str, curve)) + "]")
print("error_bar = [" + ",".join(map(str, error_bar)) + "]")
```

Then run:

```bash
cd "Code Circuits"
python RetrieveInformation.py
```

The parameters mean:

| Parameter      | Meaning                                                    |
| -------------- | ---------------------------------------------------------- |
| `filename`     | Raw or combined result file                                |
| `sweep_size`   | Number of noise points in one job; normally equal to `num` |
| `times`        | Number of repeated dictionaries used from each noise row   |
| `prepend_zero` | Adds a zero-noise point at the start of both arrays        |

The script prints:

```text
curve = [...]
error_bar = [...]
```

`curve` contains the mean logical-error probability at each noise value.

`error_bar` contains the standard error of the mean.

The current retrieval code is written for one encoded logical qubit, with logical dictionaries of the form:

```python
{"0": count_0, "1": count_1}
```

---

## Using the retrieval functions from another script

The example code at the bottom of `RetrieveInformation.py` currently runs whenever the file is imported. For library use, place it under a main guard:

```python
if __name__ == "__main__":
    fn = "test.txt"

    curve, error_bar = load_error_curve(
        filename=fn,
        sweep_size=20,
        times=10,
        prepend_zero=True,
    )

    print("curve = [" + ",".join(map(str, curve)) + "]")
    print("error_bar = [" + ",".join(map(str, error_bar)) + "]")
```

After this change, the function can be imported safely:

```python
from RetrieveInformation import load_error_curve

curve, error_bar = load_error_curve(
    filename="Enc513Full-m=10-ec=1-combined.txt",
    sweep_size=20,
    times=10,
    prepend_zero=True,
)
```

---

## Complete workflow

```text
Choose the code:
    [[5,1,3]] or [[7,1,3]]
                |
                v
Choose the model:
    Full continuous / Approximate continuous / Pauli
                |
                v
Open the matching folder
                |
                v
Edit the parameters and output filename in Tests.py
                |
                v
Run several independent parallel jobs
                |
                v
Store one raw noise sweep per job
                |
                v
Concatenate the compatible job files
                |
                v
Run RetrieveInformation.py
                |
                v
Obtain the logical-error curve and its error bars
```

---

## Reproducibility

The continuous models use NumPy to sample random angles. To reproduce the same sampled circuit realization, set a seed before the noisy circuit is built:

```python
import numpy as np

np.random.seed(1234)
```

Independent parallel jobs should use different seeds. A simple approach is to derive the seed from the job index:

```python
np.random.seed(1000 + job_index)
```

Store the seed together with the output filename or job log.

---

## Common mistakes

### Running from the wrong folder

Because the same module names are reused, running a script from another directory may load the wrong circuit code.

Use:

```bash
cd "path/to/the/chosen/model-and-code-folder"
python Tests.py
```

### Mixing incompatible raw jobs

Do not combine jobs with different values of `num`, `n_h`, `correction`, `shots`, `times`, code type, or noise model.

### Forgetting append mode

Running the same script again adds data to the existing file. Remove the file or change its name when starting a new data set.

### Using the wrong `sweep_size`

If `num=20` in `Tests.py`, use:

```python
sweep_size=20
```

Otherwise, `RetrieveInformation.py` will group the rows incorrectly or raise an error.

### Comparing raw counts instead of probabilities

The continuous and Pauli folders use different shot conventions. Compare the normalized logical-error probabilities, not the raw count totals.

### Treating the approximate model as exact

The approximate model keeps local error widths but drops multi-qubit correlations generated during propagation. Use the full continuous model when these correlations are important.

---

## Citation

This code accompanies:

> Y. El Kaderi, A. Honecker, and I. Andriyanova,
> *Continuous Noise Model for Quantum Circuits*.

Please cite the corresponding work and this repository when using the code or the generated data.
::: 

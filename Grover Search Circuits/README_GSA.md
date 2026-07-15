# Grover Search Circuits under Continuous and Pauli Noise

This folder contains Qiskit simulations of Grover's search algorithm under two noise models:

1. a continuous coherent-noise model;
2. a symmetric Pauli depolarizing channel.

The two simulation scripts have similar structures because they implement the same Grover search circuit. However, the noise is introduced in different ways, and the scripts do not use exactly the same sampling procedure.

The continuous-noise simulation was designed to be run many times in parallel. Each independent job stores one raw error-probability curve in a text file. The script `retrieveinfo.py` reads these raw curves and computes their point-by-point average.

---

## Folder structure

```text
Grover Search Circuits/
├── GSA_Continuous.py
├── GSA_Pauli.py
└── retrieveinfo.py
```

### File roles

| File                | Purpose                                                                                  |
| ------------------- | ---------------------------------------------------------------------------------------- |
| `GSA_Continuous.py` | Simulates Grover search with random coherent rotations inserted after single-qubit gates |
| `GSA_Pauli.py`      | Simulates Grover search with a symmetric Pauli channel using Qiskit Aer                  |
| `retrieveinfo.py`   | Reads and averages the raw curves produced by parallel continuous-noise runs             |

---

## Installation

The simulations require:

* Python 3.10 or later;
* NumPy;
* Matplotlib;
* Qiskit;
* Qiskit Aer.

Install the required packages with:

```bash
python -m pip install numpy matplotlib qiskit qiskit-aer
```

Move to this directory before running the scripts:

```bash
cd "Grover Search Circuits"
```

---

# Continuous coherent-noise simulation

The continuous model is implemented in:

```text
GSA_Continuous.py
```

Run it with:

```bash
python GSA_Continuous.py
```

## Circuit construction

For an `n`-qubit Grover search circuit, the script:

1. prepares the uniform superposition;
2. builds a phase oracle for the marked state;
3. builds the Grover diffuser;
4. decomposes the oracle and diffuser into lower-level gates;
5. inserts coherent errors after single-qubit gates;
6. repeats the oracle and diffuser;
7. measures all qubits;
8. computes the probability of not measuring the marked state.

The number of Grover iterations is approximated by

```python
int(np.pi / 4 * np.sqrt(2**n) - 0.5)
```

---

## Continuous error model

The coherent error applied to one qubit is decomposed as

```text
Z Rx(2 theta) Rz(-2 phi) Z
```

The angles are sampled independently from Gaussian distributions:

```text
theta ~ N(0, (k/3)^2)
phi   ~ N(0, (k/3)^2)
```

where `k` is the continuous-noise strength.

Fresh noise is inserted only after single-qubit gates. No new physical error is added directly after multi-qubit gates.

The circuit itself already contains the sampled coherent errors. It is therefore executed using an ideal `AerSimulator` without an additional Aer noise model.

---

## Main continuous-model parameters

The simulation parameters are set at the bottom of `GSA_Continuous.py`.

A typical call is:

```python
results = RunGSA_AllNoise(
    a=-3,
    b=-1,
    num=20,
    shots=200,
    n=6,
    seed=None,
)

save_table(
    [results],
    "GSA6qubits_continuous.txt",
)
```

The parameters are:

| Parameter | Meaning                                       |
| --------- | --------------------------------------------- |
| `a`       | Lower exponent of the logarithmic noise range |
| `b`       | Upper exponent of the logarithmic noise range |
| `num`     | Number of tested noise values                 |
| `shots`   | Number of 100-shot simulator batches          |
| `n`       | Number of qubits                              |
| `seed`    | Optional seed for selecting the target states |

The noise values are generated with:

```python
np.logspace(a, b, num=num)
```

For example,

```python
a = -3
b = -1
num = 20
```

produces 20 logarithmically spaced values between `10^-3` and `10^-1`.

---

## Meaning of `shots` in the continuous model

The `shots` parameter is not passed directly to Qiskit Aer as the total number of measurements.

The circuit is executed in batches of 100 shots:

```python
for _ in range(shots):
    batch_counts = QuantumSimulator(qc, shots=100)
```

Therefore,

```text
total measurements per target and noise value = 100 × shots
```

For example,

```text
shots = 200
```

corresponds to

```text
200 × 100 = 20,000 measurements
```

for each marked target and each noise value.

All these batches reuse the same noisy circuit. They improve the measurement statistics but do not generate new coherent-error angles.

A new coherent-error realization is generated when the circuit is rebuilt.

---

## Target states in the continuous script

The current `binary_strings(n)` function does not generate all `2^n` computational-basis states.

Instead, it generates `n` target states:

* one state with one zero;
* one state with two zeros;
* ...
* one state with `n` zeros.

The positions of the zeros are selected randomly.

For example, for `n=6`, the script tests six marked states rather than all 64 possible states.

To test all computational-basis targets, replace the function with:

```python
def binary_strings(n: int) -> list[str]:
    return [
        format(i, f"0{n}b")
        for i in range(2**n)
    ]
```

This provides a more complete average but increases the simulation cost.

---

# Parallel execution and raw result files

The continuous simulation requires many independent random circuit realizations. It was therefore parallelized.

Each parallel job runs the same noise sweep and generates one raw failure-probability curve. The curve is stored in a text file instead of being averaged immediately.

The saving function uses append mode:

```python
with open(filename, "a") as file:
```

This allows several independent runs to be stored in the same file.

A typical raw file has the form:

```text
Shots = 200
0.0012 0.0018 0.0025 ... 0.1840
Shots = 200
0.0010 0.0017 0.0027 ... 0.1815
Shots = 200
0.0011 0.0019 0.0026 ... 0.1862
```

Each numerical row represents one complete noise sweep from one independent run.

For `num=20`, each numerical row must contain 20 values.

---

## Recommended parallel workflow

It is safer for each parallel process to write to its own file.

For example:

```text
GSA_cont_n6_job_001.txt
GSA_cont_n6_job_002.txt
GSA_cont_n6_job_003.txt
```

After all jobs have finished, combine the files:

```bash
cat GSA_cont_n6_job_*.txt > GSA_cont_n6_all.txt
```

The combined file can then be processed with `retrieveinfo.py`.

Using separate files avoids several processes writing to the same output file at the same time.

---

## Output-file naming

Use filenames that describe the simulation settings.

A recommended format is:

```text
GSA_<model>_n=<qubits>_points=<num>_job=<id>.txt
```

For example:

```text
GSA_continuous_n=6_points=20_job=007.txt
```

Do not combine files produced with different:

* qubit numbers;
* noise ranges;
* numbers of noise points;
* shot counts;
* target-selection methods.

The retrieval script averages every valid numerical row and cannot detect whether two rows came from different experiments.

---

# Retrieving the averaged continuous result

The raw curves are processed using:

```text
retrieveinfo.py
```

Run it with:

```bash
python retrieveinfo.py
```

The script:

1. opens the selected text file;
2. skips empty lines;
3. skips lines beginning with `Shots`;
4. converts each numerical row into an array;
5. checks that each row contains the expected number of values;
6. computes the column-wise mean.

If the file contains `R` independent runs and `M` noise values, the data form an array of shape

```text
R × M
```

The final value at noise index `j` is

```text
mean[j] = average of column j over all independent runs.
```

---

## Selecting the input file

At the bottom of `retrieveinfo.py`, change the input filename:

```python
mean_values = read_mean_array(
    "GSA_cont_n6_all.txt",
    expected_length=20,
)

print(mean_values.tolist())
```

The value of `expected_length` must match `num` in `GSA_Continuous.py`.

For example, when

```python
num = 30
```

use:

```python
mean_values = read_mean_array(
    "GSA_cont_n6_all.txt",
    expected_length=30,
)
```

---

## Importing the retrieval function

The retrieval function can also be imported into another Python script:

```python
from retrieveinfo import read_mean_array

mean_curve = read_mean_array(
    "GSA_cont_n6_all.txt",
    expected_length=20,
)

print(mean_curve)
```

The current retrieval script computes only the mean. It does not compute:

* the standard deviation;
* the standard error;
* confidence intervals;
* the corresponding noise-axis values.

The noise axis can be rebuilt using the same parameters as the simulation:

```python
x = np.logspace(a, b, num=num)
```

---

# Pauli depolarizing-channel simulation

The Pauli model is implemented in:

```text
GSA_Pauli.py
```

Run it with:

```bash
python GSA_Pauli.py
```

## Pauli channel

The script applies a symmetric Pauli channel with total error probability `p`:

```text
P(X) = p/3
P(Y) = p/3
P(Z) = p/3
P(I) = 1 - p
```

The channel is added using a Qiskit Aer `NoiseModel`.

The ideal Grover circuit is constructed first. The noise model is then attached to the supported single-qubit basis gates during simulation.

---

## Main Pauli-model parameters

A typical parameter block is:

```python
results = RunGSA_AllNoise_Pauli(
    a=-6,
    b=-1,
    num=20,
    shots=100000,
    n=6,
)
```

The parameters are:

| Parameter | Meaning                                 |
| --------- | --------------------------------------- |
| `a`       | Lower exponent of the Pauli-error range |
| `b`       | Upper exponent of the Pauli-error range |
| `num`     | Number of error probabilities           |
| `shots`   | Number of Aer measurements per target   |
| `n`       | Number of qubits                        |

Unlike the continuous script, the Pauli script passes `shots` directly to Aer.

Therefore,

```text
total measurements per target and noise value = shots
```

The current Pauli script prints the result and plots the curve directly.

---

## Saving the Pauli result

To store the Pauli curve in a text file, add:

```python
np.savetxt(
    "GSA_pauli_n6.txt",
    np.asarray(results)[None, :],
)
```

This writes one row containing the complete Pauli-noise sweep.

---

# Comparing the continuous and Pauli models

Both scripts use the same performance metric:

```text
failure probability
    = 1 - probability of measuring the marked target
```

However, their default settings are not directly identical.

## Target-state averaging

The current continuous script tests `n` sampled target states.

The Pauli script tests all `2^n` computational-basis target states.

For a direct comparison, use the same target set in both simulations.

## Noise parameters

The continuous script scans `k`, while the sampled angular standard deviation is `k/3`.

The Pauli script scans the total non-identity probability `p`.

The parameters `k` and `p` should not automatically be treated as equal. Use the matching rule defined for the continuous and Pauli models before comparing the curves.

## Shot conventions

The scripts use different shot conventions:

```text
Continuous model: 100 × shots
Pauli model:      shots
```

Compare normalized probabilities rather than raw count totals.

## Independent realizations

In the continuous model, each circuit construction samples new coherent-error angles.

In the Pauli model, each shot samples the stochastic Pauli channel through Aer.

Several independent continuous runs are therefore stored and averaged to obtain a stable result.

---

# Reproducibility

The continuous simulation uses two random-number generators:

* Python's `random` module selects target states;
* NumPy samples the coherent-error angles.

The current optional `seed` controls only the Python `random` module.

To reproduce both the target selection and the coherent errors, seed both generators:

```python
import random
import numpy as np

seed = 1234

random.seed(seed)
np.random.seed(seed)
```

For parallel runs, assign a different seed to every job.

For example:

```python
seed = 1000 + job_index

random.seed(seed)
np.random.seed(seed)
```

Record the seed in the job log or output filename.

---

# Typical workflow

```text
Choose the number of qubits
            |
            v
Set the noise range and shot count
            |
            v
Run GSA_Continuous.py in several independent jobs
            |
            v
Store one raw curve per job
            |
            v
Combine the compatible raw files
            |
            v
Run retrieveinfo.py
            |
            v
Obtain the mean continuous-noise curve
            |
            v
Run GSA_Pauli.py
            |
            v
Compare the two normalized failure-probability curves
```

---

# Example continuous-model workflow

Edit `GSA_Continuous.py`:

```python
seed = 1

random.seed(seed)
np.random.seed(seed)

results = RunGSA_AllNoise(
    a=-4,
    b=-1,
    num=20,
    shots=200,
    n=6,
    seed=seed,
)

save_table(
    [results],
    "GSA_cont_n6_job_001.txt",
)
```

Run several jobs with different seeds:

```bash
python GSA_Continuous.py
```

Combine their output files:

```bash
cat GSA_cont_n6_job_*.txt > GSA_cont_n6_all.txt
```

Edit `retrieveinfo.py`:

```python
mean_values = read_mean_array(
    "GSA_cont_n6_all.txt",
    expected_length=20,
)

print(mean_values.tolist())
```

Run:

```bash
python retrieveinfo.py
```

---

# Common mistakes

## Using the wrong input filename in `retrieveinfo.py`

Change the filename at the bottom of the script to the actual combined output file.

## Using the wrong expected row length

The value of `expected_length` must equal `num` in the continuous simulation.

## Mixing incompatible runs

Do not average curves generated with different parameter sets.

## Forgetting append mode

The continuous script appends to an existing file. Delete or rename the old file before starting a different experiment.

## Assuming each 100-shot batch is a new noise realization

All batches inside `runcircuitnoisy()` reuse the same sampled noisy circuit.

## Comparing raw counts

The continuous and Pauli scripts use different shot conventions. Compare normalized failure probabilities.

## Comparing different target sets

The current continuous and Pauli scripts do not average over the same marked states by default.

---

# Citation

This code accompanies:

> Y. El Kaderi, A. Honecker, and I. Andriyanova,
> *Continuous Noise Model for Quantum Circuits*.

Please cite the corresponding work and the GitHub repository when using these scripts or their generated data.

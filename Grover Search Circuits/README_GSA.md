# Grover Search Circuits

This folder contains simulations of Grover's search algorithm under:

* continuous coherent noise;
* a symmetric Pauli depolarizing channel.

The scripts may look similar because they implement the same Grover circuit with different noise models.

## Files

```text
Grover Search Circuits/
├── GSA_Continuous.py
├── GSA_Pauli.py
└── retrieveinfo.py
```

| File                | Purpose                                                                        |
| ------------------- | ------------------------------------------------------------------------------ |
| `GSA_Continuous.py` | Simulates Grover search with coherent errors inserted after single-qubit gates |
| `GSA_Pauli.py`      | Simulates Grover search with a symmetric Pauli channel using Qiskit Aer        |
| `retrieveinfo.py`   | Reads and averages raw curves produced by repeated continuous-noise runs       |

## Installation

```bash
python -m pip install numpy matplotlib qiskit qiskit-aer
cd "Grover Search Circuits"
```

## Continuous-noise simulation

Run:

```bash
python GSA_Continuous.py
```

The coherent error applied to each noisy qubit is

```text
Z Rx(2 theta) Rz(-2 phi) Z
```

with

```text
theta, phi ~ N(0, (k/3)^2),
```

where `k` is the noise strength.

Noise is inserted after single-qubit gates in the decomposed Grover oracle and diffuser.

The main parameters are set at the bottom of the file:

```python
results = RunGSA_AllNoise(
    a=-3,
    b=-1,
    num=20,
    shots=200,
    n=6,
    seed=None,
)

save_table([results], "GSA6_continuous.txt")
```

| Parameter | Meaning                                  |
| --------- | ---------------------------------------- |
| `a`, `b`  | Exponents of the logarithmic noise range |
| `num`     | Number of tested noise values            |
| `shots`   | Number of 100-shot batches               |
| `n`       | Number of qubits                         |
| `seed`    | Optional seed for target selection       |

In this script,

```text
total measurements = 100 × shots
```

for each target and noise value.

The current target-generation function tests `n` sampled target states, not all `2^n` computational-basis states.

## Parallel runs and raw output

The continuous simulation requires many independent coherent-error realizations and was therefore parallelized.

Each job produces one complete noise-sweep curve and stores it in a text file. Append mode is used so that repeated runs can be collected:

```text
Shots = 200
0.0012 0.0018 0.0025 ... 0.1840
Shots = 200
0.0010 0.0017 0.0027 ... 0.1815
```

Each numerical row is one independent raw curve.

For parallel execution, it is safer to use one file per job:

```text
GSA_cont_n6_job_001.txt
GSA_cont_n6_job_002.txt
GSA_cont_n6_job_003.txt
```

Combine them after all jobs finish:

```bash
cat GSA_cont_n6_job_*.txt > GSA_cont_n6_all.txt
```

Only combine runs produced with the same noise range, qubit number, number of points, shot count, and target-selection method.

## Retrieving the averaged curve

Run:

```bash
python retrieveinfo.py
```

The script:

1. skips empty lines and `Shots` headers;
2. reads each numerical row;
3. checks the number of values;
4. computes the column-wise mean over all runs.

Set the input file at the bottom of `retrieveinfo.py`:

```python
mean_values = read_mean_array(
    "GSA_cont_n6_all.txt",
    expected_length=20,
)

print(mean_values.tolist())
```

`expected_length` must match `num` in `GSA_Continuous.py`.

The script returns the mean curve only. The corresponding noise axis can be rebuilt with:

```python
x = np.logspace(a, b, num=num)
```

## Pauli-noise simulation

Run:

```bash
python GSA_Pauli.py
```

The Pauli channel is

```text
P(X) = p/3
P(Y) = p/3
P(Z) = p/3
P(I) = 1 - p
```

where `p` is the total non-identity error probability.

A typical run is:

```python
results = RunGSA_AllNoise_Pauli(
    a=-6,
    b=-1,
    num=20,
    shots=100000,
    n=6,
)
```

In this script, `shots` is passed directly to Qiskit Aer.

The result is printed and plotted. It can also be saved with:

```python
np.savetxt(
    "GSA_pauli_n6.txt",
    np.asarray(results)[None, :],
)
```

## Comparing the models

Both scripts compute

```text
failure probability
    = 1 - probability of measuring the marked state.
```

Before comparing their curves, note that:

* the continuous script currently tests `n` sampled targets;
* the Pauli script tests all `2^n` targets;
* the continuous parameter `k` is not the same as the Pauli probability `p`;
* the two scripts use different shot conventions.

Use the same target set and the appropriate noise-matching rule for a direct comparison.

## Reproducibility

The continuous script uses both Python and NumPy random generators. To reproduce target selection and coherent errors, seed both:

```python
import random
import numpy as np

seed = 1234

random.seed(seed)
np.random.seed(seed)
```

Use a different seed for each parallel job.

## Citation

This code accompanies:

> Y. El Kaderi, A. Honecker, and I. Andriyanova,
> *Continuous Noise Model for Quantum Circuits* [arXiv:2604.26008](https://arxiv.org/abs/2604.26008).

Please cite the corresponding work and this repository when using the code or generated data.

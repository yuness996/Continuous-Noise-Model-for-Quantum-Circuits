# Quantum-Circuit Simulation Codes

This folder contains circuit-level simulations for the `[[5,1,3]]` and `[[7,1,3]]` quantum error-correcting codes.

Three noise models are included:

* the full continuous coherent-noise model;
* the approximate continuous-noise model;
* the Pauli depolarizing channel.

Several files have the same names because each code and noise model is stored in its own folder. Do not mix files from different folders.

---

## Folder structure

```text
Code Circuits/
├── Continuous Error Model/
│   ├── Full Model/
│   │   ├── [[5,1,3]]ContinuousError/
│   │   └── [[7,1,3]]ContinuousError/
│   │
│   └── Approximate Error Model/
│       ├── [[5,1,3]] (continuous error) (model)/
│       └── [[7,1,3]] (continuous error) (model)/
│
├── Pauli Depolarizing Channel/
│   ├── [[5,1,3]] Pauli Errors/
│   └── [[7,1,3]] Pauli Errors/
│
└── RetrieveInformation.py
```

Each simulation folder contains:

```text
Oracles.py
PostProcessing.py
Tests.py
```

Their roles are:

* `Oracles.py`: defines the encoding, logical gates, syndrome extraction, and recovery circuits;
* `PostProcessing.py`: runs the circuits and converts physical counts into logical results;
* `Tests.py`: sets the simulation parameters, runs the noise sweep, and saves the raw results.

Users normally need to edit and run `Tests.py`.

---

## Choosing a model

| Model                          | Folder                                            |
| ------------------------------ | ------------------------------------------------- |
| Full gate-level coherent noise | `Continuous Error Model/Full Model/`              |
| Faster local approximation     | `Continuous Error Model/Approximate Error Model/` |
| Symmetric Pauli channel        | `Pauli Depolarizing Channel/`                     |

The full continuous model inserts random coherent errors directly after noisy gates.

The approximate model replaces the accumulated gate errors by one final effective error on each qubit. It is faster, but it does not keep all multi-qubit error correlations.

The Pauli model applies

```text
P(X) = P(Y) = P(Z) = p/3,
P(I) = 1 - p.
```

---

## Installation

Install the required packages with:

```bash
python -m pip install numpy sympy qiskit qiskit-aer
```

`SymPy` is needed for the approximate model.

---

## Running a simulation

Move to the folder of the selected code and model, then run:

```bash
python Tests.py
```

For example:

```bash
cd "Code Circuits/Continuous Error Model/Full Model/[[5,1,3]]ContinuousError"
python Tests.py
```

Run the script from its own folder so that it imports the correct local versions of `Oracles.py` and `PostProcessing.py`.

The main parameters are set near the end of `Tests.py`.

A typical call is:

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

The main parameters are:

| Parameter    | Meaning                                      |
| ------------ | -------------------------------------------- |
| `a`, `b`     | Exponents of the logarithmic noise range     |
| `num`        | Number of tested noise values                |
| `n_h`        | Number of logical Hadamard gates             |
| `correction` | Enables syndrome extraction and recovery     |
| `shots`      | Number of measurement samples                |
| `times`      | Number of repeated runs for each noise value |

Check the selected `Tests.py`, since some function names and parameters differ between folders.

---

## Parallel simulations and raw files

The simulations require many noisy circuit runs. They were therefore executed as independent parallel jobs.

Each job performs a full sweep over the noise values and saves the raw logical counts in a text file. The files are stored instead of directly saving the final curve so that results from several jobs can later be combined.

The scripts use append mode:

```python
with open(filename, "a") as file:
```

Running the same script again adds new results to the existing file.

Use a new output name when starting a different experiment. Do not mix results produced with different codes, noise models, shot counts, circuit depths, or error-correction settings.

For parallel runs, it is safer to use one file per job:

```text
513_full_m10_ec1_job001.txt
513_full_m10_ec1_job002.txt
513_full_m10_ec1_job003.txt
```

After all jobs finish, combine them:

```bash
cat 513_full_m10_ec1_job*.txt > 513_full_m10_ec1_all.txt
```

---

## Retrieving the final results

Use:

```text
RetrieveInformation.py
```

This script reads the raw files produced by the parallel jobs. It combines matching noise points, computes the logical error probability, and returns the mean curve and its error bars.

Set the input file near the end of the script:

```python
filename = "513_full_m10_ec1_all.txt"

curve, error_bar = load_error_curve(
    filename=filename,
    sweep_size=20,
    times=10,
    prepend_zero=True,
)
```

Then run:

```bash
cd "Code Circuits"
python RetrieveInformation.py
```

The parameter `sweep_size` must match `num` in `Tests.py`.

The script prints:

```text
curve = [...]
error_bar = [...]
```

These arrays can then be used for plotting.

---

## Important notes

The continuous and Pauli scripts do not always use `shots` in the same way. In the continuous folders, one value of `shots` may correspond to several batches of Aer shots. In the Pauli folders, it is usually passed directly to Aer.

Compare normalized logical error probabilities, not raw count totals.

The continuous model also samples random errors. Set a NumPy seed when a reproducible run is needed:

```python
import numpy as np

np.random.seed(1234)
```

Use a different seed for each independent parallel job.

---

## Citation

This code accompanies:

> Y. El Kaderi, A. Honecker, and I. Andriyanova,
> *Continuous Noise Model for Quantum Circuits* [arXiv:2604.26008](https://arxiv.org/abs/2604.26008).

Please cite the corresponding work and this repository when using the code or the generated data.
::: 

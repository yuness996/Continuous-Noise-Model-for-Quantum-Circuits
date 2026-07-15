# Grover Search under Continuous and Pauli Noise

This folder contains two Qiskit simulations of Grover's search algorithm. Both scripts compute the probability of failing to measure the marked computational-basis state, but they use different physical noise models.

## Files

### `GSA_Continuous.py`

Simulates Grover search with random coherent one-qubit rotations. After each noisy one-qubit gate, the circuit applies

\[
U_{\mathrm{err}}(\theta,\phi)
= ZR_x(2\theta)R_z(-2\phi)Z,
\]

where the two angles are sampled independently from

\[
\theta,\phi \sim \mathcal{N}\left(0,\frac{k^2}{9}\right).
\]

The phase oracle and diffuser are decomposed before the noise gates are inserted. No new physical noise is attached directly to multi-qubit gates.

For each noise strength, the script:

1. selects one representative target for each possible number of zero bits;
2. builds a noisy Grover circuit for each target;
3. runs the same sampled circuit in several batches of 100 shots;
4. merges the counts; and
5. averages the failure probability over the selected targets.

A key point is that the random angles are sampled when the circuit is built. All shot batches for that circuit therefore use the same coherent-error realization. To average over independent noise realizations, the circuit must be rebuilt inside the batch loop.

### `GSA_Pauli.py`

Simulates Grover search with a symmetric stochastic Pauli channel. After each supported one-qubit basis gate, Aer applies

\[
\mathcal{E}(\rho)
= (1-p)\rho
+ \frac{p}{3}X\rho X
+ \frac{p}{3}Y\rho Y
+ \frac{p}{3}Z\rho Z.
\]

For each value of `p`, the script runs the circuit for all \(2^n\) possible marked states and reports the average failure probability.

## Performance metric

For a marked state \(\lvert\omega\rangle\), the measured failure probability is

\[
p_{\mathrm{err}}^{\mathrm{GSA}}
= 1-p_{\lvert\omega\rangle},
\]

where \(p_{\lvert\omega\rangle}\) is the observed probability of the target bitstring. Equivalently,

\[
p_{\mathrm{err}}^{\mathrm{GSA}}
= 1-\frac{N_{\omega}}{N_{\mathrm{shots}}}.
\]

The final value at each noise strength is averaged over the tested target states.

## Requirements

- Python 3.9 or later
- NumPy
- Matplotlib
- Qiskit
- Qiskit Aer

Install the dependencies with

```bash
python -m pip install numpy matplotlib qiskit qiskit-aer
```

## Running the simulations

Run the continuous-noise model with

```bash
python GSA_Continuous.py
```

The default configuration uses six qubits, 20 logarithmically spaced noise strengths from \(10^{-3}\) to \(10^{-1}\), and 200 batches of 100 shots. Results are appended to `GSA6qubits_ed2.txt`.

Run the Pauli model with

```bash
python GSA_Pauli.py
```

The default configuration uses six qubits, all 64 target states, 20 logarithmically spaced error probabilities from \(10^{-6}\) to \(10^{-1}\), and 100,000 shots per target. The result array is printed and plotted.

## Main parameters

In both scripts:

- `n_qubits` sets the search-space size to \(2^{n}\).
- `exponent_min` and `exponent_max` set the logarithmic noise range.
- `num_points` sets the number of noise values.
- `shots` or `batches` controls the sampling cost.

The number of Grover iterations is

\[
r=\left\lfloor \frac{\pi}{4}\sqrt{2^n}-\frac{1}{2}\right\rfloor.
\]

## Model-comparison notes

The parameters `k` and `p` are not directly the same quantity. In the continuous model, `k/3` is the standard deviation of each sampled angle. In the Pauli model, `p` is the total probability of a non-identity Pauli error after a noisy one-qubit gate. A physical comparison requires a stated mapping, such as matching the average channel infidelity or another common error metric.

The target averaging also differs by default. The Pauli script uses every computational-basis target, while the continuous script uses only `n_qubits` sampled targets. Use the same target set in both scripts when making a direct numerical comparison.

## Reproducibility

The continuous script accepts a random seed in `run_grover_noise_sweep`. It seeds both Python's `random` module and NumPy's random generator. The Pauli simulation remains subject to finite-shot fluctuations unless an Aer simulator seed is also supplied.

## Repository

The complete project is available at:

`https://github.com/yuness996/Continuous-Noise-Model-for-Quantum-Circuits`

## Processing repeated simulation results

### `retrieveinfo_commented.py`

This helper script reads a text file containing several runs of the same GSA
noise sweep and computes the mean result at each noise point. It is intended
for output files produced by the continuous-noise script when several result
arrays are appended to the same file.

The function

```python
read_mean_array(filename, expected_length=20)
```

performs the following steps:

1. opens the selected text file;
2. ignores blank lines and metadata lines beginning with `Shots`;
3. converts each remaining line into floating-point values;
4. checks that every numerical row contains the expected number of points; and
5. averages the rows column by column with `numpy.mean(..., axis=0)`.

If the file contains result vectors

\[
\mathbf{p}^{(1)},\mathbf{p}^{(2)},\ldots,\mathbf{p}^{(R)},
\]

the returned array is

\[
\overline{p}_j=\frac{1}{R}\sum_{r=1}^{R}p_j^{(r)},
\]

where \(j\) labels the noise-strength point and \(R\) is the number of stored
runs. Thus, the script averages repeated estimates at the same physical noise
strength; it does not average across different noise strengths.

The expected input format is

```text
Shots = 1000
0.001 0.002 ... 0.020
Shots = 1000
0.0011 0.0019 ... 0.0198
```

By default, each numerical line must contain 20 values, matching the default
20-point logarithmic sweep used by the GSA scripts. A different sweep length
can be supplied explicitly:

```python
mean_values = read_mean_array("results.txt", expected_length=30)
```

Run the default example with

```bash
python retrieveinfo_commented.py
```

The script reads `GSA7qubits_ed2.txt` from the current working directory and
prints the column-wise mean as a standard Python list. The input filename can
be changed in `main()`, or the function can be imported into another analysis
script.

The script reports a clear error when the input file is missing, contains
non-numeric data, has rows of inconsistent length, or contains no numerical
results.

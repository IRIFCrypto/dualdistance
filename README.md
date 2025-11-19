Bounding the Dual Distance of Sparse Matrices
=============================================

This is the artifact of the paper "Fast Pseudorandom Correlation Functions from Sparse LPN" by
Lennart Braun, Geoffroy Couteau, Kelsey Melissaris, Mahshid Riahinia, and Elahe Sadeghi
which will be published at Asiacrypt 2025.
The full version is available on [IACR ePrint](https://eprint.iacr.org/2025/1644).

The artifact implements the method described in Section 4.3 to compute a lower bound $D$ on the dual
distance of random $k$-sparse regular $m \times n$ matrices $\mathbf{A}$ over $\mathbb{F}_2$ where
each row is the concatenation of $k$ random unit vectors of length $N := n/k$.
The bound $D$ holds except with probability at most $2^{-\rho}$, where $\rho$ is a statistical
security parameter.
It has been used to obtain the values given in Table 2.


Mistakes in the Conference Version
---

A previous version of the script, which was used to compute the values of $D$ reported in the
conference version and the current ePrint version, contained two mistakes:

- The reported value for $D$ was too small: it was given as the largest $w$ such that there is no
  attack vector of Hamming weight $\leq w2$ instead of being a strict upper bound on the Hamming
  weight $2w$ itself.
- The computation of $p_{w,w}$ was incorrect, resulting in slightly too large values.

Combined the two errors resulted in values for $D$ that were too small by almost exactly a factor 2.
Since $D$ is a lower bound, the reported values were *not* incorrect, but not as good as they should
have been.

While the ePrint version has not been corrected yet, we will provide an new version of the paper
with the artifact submission where Table 2 and the resulting PCF parameters have been updated.
It can be verified against the ePrint version that the new values for $D$ are twice as large than
before.


Dependencies
------------

- Python 3
- [`gmpy2`](https://github.com/aleaxit/gmpy)
- [`tqdm`](https://github.com/tqdm/tqdm) (optional, if installed shows a progress bar)


Running
-------

Parameters are passed as `key=value` pairs on the command line:
- `statsecs`: comma-separated integers, statistical security parameters
- `log_n`: integer, base-2 logarithm of the number of columns $n$
- `k`: integer, sparsity of the matrix
- `s`: integer, expansion factor, i.e., ratio $s = m/n$
- `cutoff`: integer, compute the union bound only up to $w < N/cutoff$
- `verbose`: integer, optional, output extra information to stderr if non-zero
- `json`: integer, optional, output results in JSON if non-zero


Output
------

For each statistical security parameter, the script will output values for
- `D`: $D$
- `delta`: $\delta = D/n$
- `union_bound`: the last value of the union bound before passing the $2^{-\rho}$ threshold

If no non-trivial bound was found, all three values are zero.
If the `cutoff` parameter results in the script stopping before passing the threshold, a value
`'stopped_early': True` is included in the results.


Examples
--------

The following examples show how to recompute single entries in Table 2.


### With Verbose Output

```
python compute_D.py statsecs=40,64,80 log_n=16 k=10 s=32 cutoff=2 verbose=1
[+] Computing transition probabilities ...
100%|███████████████████████████████████████████████████████████████████████| 3276/3276 [00:00<00:00, 100238.84it/s]
[+] Main loop over w ...
 75%|███████████████████████████████████████████████████████▍                  | 2451/3275 [00:05<00:01, 448.84it/s][-] sum_(w=1)^2478 q_w > 2^-80, stopping
[-] sum_(w=1)^2480 q_w > 2^-64, stopping
[-] sum_(w=1)^2482 q_w > 2^-40, stopping
 76%|████████████████████████████████████████████████████████                  | 2480/3275 [00:05<00:01, 440.18it/s]
{'N': 6553,
 'cutoff': 2,
 'k': 10,
 'log_n': 16,
 'm': 2097152,
 'n': 65536,
 'results': {40: {'D': 4964,
                  'delta': 0.07574462890625,
                  'union_bound': 4.85915345375458e-15},
             64: {'D': 4960,
                  'delta': 0.07568359375,
                  'union_bound': 8.995818655595886e-21},
             80: {'D': 4956,
                  'delta': 0.07562255859375,
                  'union_bound': 3.18524087717089e-26}},
 's': 32}
```

The resulting values 
$D = 4964$ for $\rho = 40$,
$D = 4960$ for $\rho = 64$,
and $\delta \approx 0.0757$
match the corresponding entries of Table 2 for $n = 2^{16}$ and $k = 10$.


### With JSON Output

(pretty-printed with [`jq`](https://jqlang.github.io/jq/))
```
python compute_D.py statsecs=40,64,80 log_n=17 k=8 s=32 cutoff=2 json=1 | jq
{
  "s": 32,
  "N": 16384,
  "cutoff": 2,
  "log_n": 17,
  "n": 131072,
  "k": 8,
  "m": 4194304,
  "results": {
    "80": {
      "D": 0,
      "union_bound": 0,
      "delta": 0.0
    },
    "64": {
      "D": 7812,
      "union_bound": 2.3149541631846313E-21,
      "delta": 0.059600830078125
    },
    "40": {
      "D": 7818,
      "union_bound": 5.18776107516655E-15,
      "delta": 0.0596466064453125
    }
  }
}
```

The resulting values
$D = 7818$ for $\rho = 40$,
$D = 7812$ for $\rho = 64$,
and $\delta \approx 0.0596$
match the corresponding entries of Table 2 for $n = 2^{17}$ and $k = 8$.

In this case no non-trivial bound $D$ was found for $\rho = 80$.

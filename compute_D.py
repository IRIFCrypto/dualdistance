from pprint import pprint
import json
import gmpy2
from gmpy2 import mpfr, bincoef, log2
import sys

# Import tqdm for progress bars if available.
try:
    from tqdm import tqdm as tqdm_func
except:
    tqdm_func = lambda x: x

# Set high precision.
gmpy2.get_context().precision = 100

def compute_transition_probs(N, cutoff, verbose=False):
    """
    Precompute transition probabilities c_(u,-1), c_(u,0), c_(u,+1) for given N and u \\in [0,cutoff+1]:
    - c_(u,-1) = bincoef(N-2, 2*(u-1)) / bincoef(N, 2*(u-1))
    - c_(u,0) = 2 * bincoef(N-2, 2*u-1) / bincoef(N, 2*u)
    - c_(u,+1) = bincoef(N-2, 2*u) / bincoef(N, 2*(u+1))
    """

    # Enable progress bars if verbose is True.
    tqdm = tqdm_func if verbose else (lambda x: x)

    # Prepare empty lists (indexed by u).
    c_um1 = [mpfr(0)] * (cutoff + 2)
    c_u0 = [mpfr(0)] * (cutoff + 2)
    c_up1 = [mpfr(0)] * (cutoff + 2)

    # Fill with the initial values for u = 0 and u = 1.
    # (NB: We set c_(0,-1) = c_(0,0) = 0 if the lower part of the binomial coefficient is negative.)
    c_um1[1] = mpfr(1)
    c_u0[1] = mpfr(2 * bincoef(N - 2, 2 * 1 - 1)) / mpfr(bincoef(N, 2 * 1))
    c_up1[0] = mpfr(bincoef(N - 2, 2 * 0)) / mpfr(bincoef(N, 2 * (0 + 1)))
    c_up1[1] = mpfr(bincoef(N - 2, 2 * 1)) / mpfr(bincoef(N, 2 * (1 + 1)))

    # In the following, we iterate over u \in [2, cutoff+1] and compute the three transition
    # probabililities for each u. Each of these probabilities is a fraction of binomial coefficients
    # (multiplied by 2 in the case of c_(u,0)), and the coefficients we need follow two sequences:
    # - In the numerators: bincoef(N-2, t) for t >= 1
    # - In the denominators: bincoef(N-2, 2*t) for t >= 1

    # To avoid computing every binomial coefficient from scratch, we use two helper functions
    # to iteratively compute the next value in two sequences of binomial coefficients.

    def next_bc1(t, prev):
        """
        Compute bincoef(N-2, t) given prev = bincoef(N, t-1).
        """
        return prev * (N - 1 - t) // t

    # Constants used in next_bc2.
    bc_t1 = N**2 + 3*N + 2
    bc_t2 = 4*N + 6
    def next_bc2(t, prev):
        """
        compute bincoef(N, 2*t) given prev = bincoef(N, 2*(t-1))
        """
        return prev * (bc_t1 + 4*t**2 - t * bc_t2) // (4*t**2 - 2*t)

    # For a given u, we need three successive values of each sequence:
    # - In the numerators: t = [2*u-2, 2*u-1, 2*u]
    # - In the denominators: t = [u-1, u, u+1]    (or equivalently 2*t = [2*u-2, 2*u, 2*u+2])
    # We store the these in the following lists, starting with the initial values for u=2:
    bcs_num = [bincoef(N-2, 2), bincoef(N-2, 3), bincoef(N-2, 4)]
    bcs_denom = [bincoef(N, 2), bincoef(N, 4), bincoef(N, 6)]

    if verbose:
        print("[+] Computing transition probabilities ...", file=sys.stderr)
    for u in tqdm(range(2, cutoff + 2)):
        # Compute c_(u,-1), c_(u,0), c_(u,+1):
        c_um1[u] = mpfr(bcs_num[0]) / mpfr(bcs_denom[0])
        c_u0[u] = mpfr(2 * bcs_num[1]) / mpfr(bcs_denom[1])
        c_up1[u] = mpfr(bcs_num[2]) / mpfr(bcs_denom[2])

        # Shift the two sequences of binomial coefficients for the next values of u.
        # - Shift the numerator sequence by two positions:
        bcs_num[0] = bcs_num[2]
        bcs_num[1] = next_bc1(2 * u + 1, bcs_num[0])
        bcs_num[2] = next_bc1(2 * u + 2, bcs_num[1])
        # - Shift the denominator  sequence by one position:
        bcs_denom = [bcs_denom[1], bcs_denom[2], next_bc2(u + 2, bcs_denom[2])]

    return c_um1, c_u0, c_up1


def compute_D(N, k, m, cutoff, statsecs, verbose=False):
    """
    Compute a lower bound D on the dual distance of a random k-sparse regular m x (kN) matrix A with
    statistical security parameter statsec.

    - N: blocksize
    - k: sparsity
    - m: number of rows
    - cutoff: number of terms to consider
    - statsec: statistical security parameter

    We find D as the largest value such that a union bound over w \\in [1,floor((D-1)/2)] of the
    probability that there exists a weight-2w vector cancelling A stays below 2^-statsec.
    """

    # Enable progress bars if verbose is True.
    tqdm = tqdm_func if verbose else (lambda x: x)

    assert len(set(statsecs)) == len(statsecs), "statsecs should be distinct"
    statsecs = sorted(statsecs, reverse=True)

    # Precompute some constants 1/N and 1-1/N.
    N_inv = mpfr(1) / mpfr(N)
    N_inv_sq = N_inv * N_inv
    one_m_N_inv = 1 - N_inv

    # Precompute the transition probabilities c_(u,-1), c_(u,0), c_(u,+1) w.r.t. N
    # for u \in [0,cutoff+1].
    c_um1, c_u0, c_up1 = compute_transition_probs(N, cutoff, verbose=verbose)

    # We want to stop as soon as the union bound exceeds the following probability:
    ub_limits = [mpfr(2) ** -s for s in statsecs]
    # Index of the first bound to check.
    statsec_i = 0

    # Recall from the paper that the value p_(w,u) is defined as the probability that the xor of 2*w
    # random unit vectors of length N := n/k has Hamming weight 2*u, and we are interested in the
    # values p_(w,0) up to some cutoff point w_max = `cutoff`.
    w_max = cutoff

    # Consider a (virtual) table with entry p_(w,u) in row w and column u, and note that the
    # recurrence can be used to compute the entry p_(w,u) from p_(w-1,u-1), p_(w-1,u), p_(w-1,u+1).
    # We do not need to compute the full table since many entries are either 0 or not necessary to
    # compute the p_(w,0) up to w = cutoff.
    # - Some entries are (considered to be) 0:
    #   - p_(w,u) = 0 for u < 0
    #   - p_(w,u) = 0 for u > w
    # - Some entries are unused:
    #   - p_(w_max,0) depends on the first i entries p_(w_max-i,0), ..., p_(w_max-i,i) of the ith
    #     previous row.
    #     => We do not need to compute p_(w_max-i,i+j) for i,j > 0.
    # Hence, we end up with a triangle shape of entries that we actually need to compute:
    #   - p_(w,u) is relevant if 0 <= w <= w_max and 0 <= u <= min(w, w_max-w)
    # This is maximized when w = w_max/2 and our table has the final shape of w_max x (w_max/2).

    # - u = 0
    # - 1 <= u <= min(w-1, w_max-w+1)
    # - u = w

    # Store the values p_(w,u) for the current iteration and  p_(w-1,u) of the previous iteration.
    ps_prev = [mpfr(0)] * (w_max // 2 + 2)
    ps = [mpfr(0)] * (w_max // 2 + 2)

    # Initial conditions
    ps[0] = N_inv
    ps[1] = (mpfr(N) - 1) / mpfr(N)

    # To compute the term q_w, we need bincoef(m, 2w). To avoid computing all the binomial
    # coefficients from scratch, we compute them iteratively with the following function.

    # Constants used in next_bc2.
    bc_t1 = m**2 + 3*m + 2
    bc_t2 = 4*m + 6
    def next_bc2(t, prev):
        """
        Compute bincoef(m, 2*t) given bincoef(m, 2*(t-1)).
        """
        return prev * (bc_t1 + 4*t**2 - t * bc_t2) // (4*t**2 - 2*t)

    # Current binomcial coefficient in the sequence, starting with w=1:
    bc_w = bincoef(m, 2) # bincoef(m, 2*w)

    # Sum of the q_i for i \in [1,w]:
    union_bound = bc_w * ps[0] ** k

    results = {}

    # Check if we are already exceeding the probability bound. In this case, we cannot find a
    # meaningful lower bound on the dual distance.
    while union_bound > ub_limits[statsec_i]:
        if verbose:
            print(f'[-] q_1 > 2^-{statsecs[statsec_i]}, no non-trivial lower bound D found', file=sys.stderr)
        results[statsecs[statsec_i]] = {
            'D': 0,
            'union_bound': 0,
        }
        statsec_i += 1
        if statsec_i >= len(statsecs):
            return results

    # For every w as long as ub_limits is not passed, we get a bound D = 2 * w + 2 because D should
    # be even and strictly larger than 2*w.
    D = 2 * 1 + 2

    # Fill in the table using the recurrence
    if verbose:
        print("[+] Main loop over w ...", file=sys.stderr)
    for w in tqdm(range(2, w_max+1)):
        # Swap the lists such that `ps_prev` contains the values from the last iterations and we can
        # overwrite `ps`.
        ps_prev, ps = ps, ps_prev

        # Use the recurrence to compute the relevant values p_(w,u) based on the previously computed
        # p_(w-1,u).

        # Fill ps with p_(w,u) for u \in [0, max(w,w_max-w)].

        # Every entry p_(w,u) depends on the three values p_(w-1,u-1), p_(w-1,u), p_(w-1,u+1), and
        # we handle the edge cases separately.

        # I. Case u = 0:
        ps[0] = N_inv * ps_prev[0] + one_m_N_inv * (
            ps_prev[1] * c_up1[0]
            # NB: the remaining terms of this sum are 0 since c_(0,-1) = c_(0,0) = 0.
        )

        # Compute the binomial coefficient
        bc_w = next_bc2(w, bc_w) # bincoef(m, 2 * w)
        # and the next term q_w of the sum:
        q_w = bc_w * ps[0] ** k

        # Check if the new term would cause us to go over the bound
        while union_bound + q_w > ub_limits[statsec_i]:
            if verbose:
                print(f'[-] sum_(w=1)^{w} q_w > 2^-{statsecs[statsec_i]}, stopping', file=sys.stderr)
            results[statsecs[statsec_i]] = {
                'D': D,
                'union_bound': float(union_bound),
            }
            statsec_i += 1
            if statsec_i >= len(statsecs):
                return results

        D = 2 * w + 2
        union_bound += q_w

        # II. Case u \in [1, min(w-1,w_max-w))]:
        for u in range(1, min(w-1,w_max-w) + 1):
            term1 = N_inv * ps_prev[u]
            term2 = one_m_N_inv * (
                ps_prev[u - 1] * c_um1[u] +
                ps_prev[u] * c_u0[u] +
                ps_prev[u + 1] * c_up1[u]
            )
            ps[u] = term1 + term2
            # XXX: maybe Eq (7) is faster

        # III. Case u = w: product of (N - i)/N from i = 1 to 2w-1
        if w <= w_max-w:
            ps[w] = ps_prev[w-1] * mpfr(N - w + 2) * mpfr(N - w + 1) * N_inv_sq

    print(f'[-] stopping early at {w=}', file=sys.stderr)
    for i in range(statsec_i, len(statsecs)):
        results[statsecs[i]] = {
            'D': D,
            'union_bound': float(union_bound),
            'stopped_early': True,
        }
    return results


def main(argv):
    kwargs = {k: v for k, v in map(lambda x: x.split('='), argv[1:])}

    if 'statsecs' not in kwargs:
        print("'statsecs' argument missing", file=sys.stderr)
        exit(1)
    try:
        statsecs = sorted(int(s) for s in kwargs['statsecs'].split(','))
        assert len(set(statsecs)) == len(statsecs)
    except:
        print("could not parse 'statsecs', it should be a comma separated list of distinct integers"
              "(e.g. 'statsecs=40,64,80')", file=sys.stderr)
        exit(1)
    if 'log_n' not in kwargs:
        print("'log_n' (dimension) argument missing", file=sys.stderr)
        exit(1)
    n = 2 ** int(kwargs['log_n'])
    if 'k' not in kwargs:
        print("'k' (sparsity) argument missing", file=sys.stderr)
        exit(1)
    k = int(kwargs['k'])
    N = n // k
    if 's' not in kwargs:
        print("'s' (stretch) argument missing", file=sys.stderr)
        exit(1)
    s = int(kwargs['s'])
    m = n * s
    if 'cutoff' not in kwargs:
        print("'cutoff' argument missing", file=sys.stderr)
        exit(1)
    cutoff = int(kwargs['cutoff'])

    verbose = 'verbose' in kwargs and int(kwargs['verbose']) != 0
    print_json = 'json' in kwargs and int(kwargs['json']) != 0

    results = compute_D(N, k, m, N // cutoff, statsecs, verbose=verbose)

    data = {}
    data['s'] = s
    data['N'] = N
    data['cutoff'] = cutoff
    data['log_n'] = int(kwargs['log_n'])
    data['n'] = n
    data['k'] = k
    data['m'] = m
    data['results'] = results
    for s in statsecs:
        data['results'][s]['delta'] = data['results'][s]['D'] / n
    if print_json:
        print(json.dumps(data))
    else:
        pprint(data)



if __name__ == "__main__":
    main(sys.argv)

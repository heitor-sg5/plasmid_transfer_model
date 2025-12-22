import numpy as np

base_params = {
    'K': 1000,       # carrying capacity
    'r': 1.0,        # baseline growth rate
    'mu': 0.1,       # baseline mortality
    'c': 0.05,       # plasmid cost
    'beta': 0.01,    # conjugation rate
    'delta': 0.01,   # plasmid loss
    's': 0.0,        # selective pressure
}

TMAX = 150
y0 = [0.9, 0.1]

def gillespie_ssa(p, TMAX, y0):
    F = int(y0[0] * p['K'])
    P = int(y0[1] * p['K'])
    t = 0.0

    times = [t]
    Fs = [F]
    Ps = [P]

    while t < TMAX and (F + P) > 0:
        N = F + P

        rF = p['r'] * (1 - p['s'])
        rP = p['r'] * (1 - p['c'])

        wF = max(0.0, rF * (1 - N / p['K']))
        wP = max(0.0, rP * (1 - N / p['K']))

        a = np.array([
            wF * F,                     # F birth
            wP * P,                     # P birth
            p['mu'] * F,                # F death
            p['mu'] * P,                # P death
            p['delta'] * P,             # plasmid loss
            p['beta'] * F * P / p['K']  # conjugation
        ])

        a0 = a.sum()
        if a0 <= 0:
            break

        tau = np.random.exponential(1 / a0)
        t += tau

        r = np.random.rand() * a0
        cumulative = np.cumsum(a)

        if r < cumulative[0]:
            F += 1
        elif r < cumulative[1]:
            P += 1
        elif r < cumulative[2]:
            F -= 1
        elif r < cumulative[3]:
            P -= 1
        elif r < cumulative[4]:
            P -= 1
            F += 1
        else:
            F -= 1
            P += 1

        F = max(F, 0)
        P = max(P, 0)

        times.append(t)
        Fs.append(F)
        Ps.append(P)

    return np.array(times), np.array(Fs), np.array(Ps)

def run_multiple_ssa(s_values=[0.0, 0.2, 0.6], n_runs=50):
    results = {}
    t_grid = np.linspace(0, TMAX, 500)

    for s in s_values:
        Fs_runs = []
        Ps_runs = []

        for _ in range(n_runs):
            p = base_params.copy()
            p['s'] = s
            t, F, P = gillespie_ssa(p, TMAX, y0)

            F_interp = np.interp(t_grid, t, F)
            P_interp = np.interp(t_grid, t, P)

            Fs_runs.append(F_interp / p['K'])
            Ps_runs.append(P_interp / p['K'])

        results[s] = {
            't': t_grid,
            'F': np.array(Fs_runs),
            'P': np.array(Ps_runs)
        }

    return results

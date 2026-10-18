"""Synthetic cooling curves for testing the analyser before real data exists.

Shape: Newton cooling towards the bath, supercooling below the freezing point,
a sharp rise when ice nucleates (recalescence), then a slowly falling plateau
(the liquid freeze-concentrates, so a solution's plateau slopes downwards),
then cooling towards the bath once mostly frozen. Probe noise and the
0.0625 C quantisation of a 12-bit DS18B20 are included.
"""

import numpy as np


def simulate(T_freeze_C, T_start_C=20.0, T_bath_C=-15.0, supercool_C=3.0,
             tau_s=180.0, dt_s=2.0, plateau_s=600.0, plateau_drop_C=1.0,
             noise_C=0.03, quantise_C=0.0625, seed=0):
    rng = np.random.default_rng(seed)
    t, T = [0.0], [T_start_C]
    target = T_freeze_C - supercool_C
    # 1. cool until supercooled
    while T[-1] > target:
        T.append(T[-1] + (T_bath_C - T[-1]) * dt_s / tau_s)
        t.append(t[-1] + dt_s)
    # 2. recalescence: jump to just below T_freeze in ~10 s
    for k in range(1, 6):
        T.append(target + (T_freeze_C - 0.05 - target) * k / 5)
        t.append(t[-1] + dt_s)
    # 3. sloping plateau
    n_p = int(plateau_s / dt_s)
    start = T[-1]
    for k in range(1, n_p + 1):
        T.append(start - plateau_drop_C * (k / n_p) ** 1.5)
        t.append(t[-1] + dt_s)
    # 4. cool to bath
    for _ in range(int(3 * tau_s / dt_s)):
        T.append(T[-1] + (T_bath_C - T[-1]) * dt_s / tau_s)
        t.append(t[-1] + dt_s)
    T = np.array(T) + rng.normal(0, noise_C, len(T))
    if quantise_C:
        T = np.round(T / quantise_C) * quantise_C
    return np.array(t), T

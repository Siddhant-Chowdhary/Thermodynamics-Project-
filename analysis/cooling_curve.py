"""Extract the freezing point from a logged cooling curve.

Method (standard for solutions):
  1. Smooth the trace lightly.
  2. Find the supercooling minimum: the lowest point that is followed by a rise
     of at least `min_jump_C` (ice nucleation releasing latent heat).
  3. The freezing point is the MAXIMUM of the smoothed trace in a window after
     that minimum. For a solution the plateau slopes downwards as the liquid
     freeze-concentrates, so its highest point is the best estimate of the
     equilibrium freezing temperature of the starting composition.
  4. Uncertainty combines probe noise (scatter of raw points near the maximum)
     and a calibration offset uncertainty supplied by the user.
"""

import numpy as np


def _smooth(y, k):
    if k <= 1:
        return y.copy()
    pad = k // 2
    ypad = np.pad(y, pad, mode="edge")
    return np.convolve(ypad, np.ones(k) / k, mode="valid")


def find_freezing_point(t, T, smooth_pts=5, min_jump_C=0.3, window_s=240.0,
                        calib_offset_C=0.0, calib_unc_C=0.1):
    t = np.asarray(t, float)
    T = np.asarray(T, float) - calib_offset_C
    Ts = _smooth(T, smooth_pts)

    # supercooling minimum followed by a clear rise
    i_min = None
    for i in range(1, len(Ts) - 1):
        ahead = Ts[i:i + max(3, int(60 / max(np.median(np.diff(t)), 1e-9)))]
        if ahead.max() - Ts[i] >= min_jump_C and Ts[i] <= Ts[i - 1]:
            # take the local minimum around i
            j_end = i + int(np.argmax(ahead))
            i_min = i + int(np.argmin(Ts[i:j_end + 1]))
            break
    if i_min is None:
        raise ValueError("No recalescence found: no rise of at least "
                         f"{min_jump_C} C after a minimum. Check the trace.")

    in_win = (t >= t[i_min]) & (t <= t[i_min] + window_s)
    idx = np.where(in_win)[0]
    i_max = idx[int(np.argmax(Ts[idx]))]
    Tf = float(Ts[i_max])

    near = np.abs(t - t[i_max]) <= 20.0
    noise = float(np.std(T[near], ddof=1)) if near.sum() > 2 else 0.05
    unc = float(np.sqrt(noise ** 2 + calib_unc_C ** 2))

    return {"T_freeze_C": Tf, "uncertainty_C": unc,
            "T_supercool_C": float(Ts[i_min]),
            "supercooling_C": Tf - float(Ts[i_min]),
            "t_nucleation_s": float(t[i_min]), "t_freeze_s": float(t[i_max])}


def load_csv(path, column=1):
    """Load a logger CSV: time_s, T1[, T2, T3]. Returns (t, T) for one probe."""
    data = np.genfromtxt(path, delimiter=",", names=True)
    names = data.dtype.names
    return data[names[0]], data[names[column]]

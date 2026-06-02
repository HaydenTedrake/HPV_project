"""
Read NANCY_pm_errs.csv, linearly interpolate PM errors vs Gaia G, plot for G ≤ 20.5,
and print one example magnitude.
"""

from pathlib import Path
import csv
import matplotlib.pyplot as plt
import numpy as np

CSV_PATH = Path(__file__).resolve().parent / "NANCY_pm_errs.csv"
G_MAX = 20.5
G_PLUG = 19.9744
N_CURVE = 800

def _parse_float(s: str) -> float:
    t = s.strip().lower()
    if t in ("inf", "+inf", "infinity"):
        return np.inf
    return float(t)

def load_csv(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    g_list, gaia_list, nancy_list = [], [], []
    with path.open(newline="", encoding="utf-8") as f:
        r = csv.reader(f)
        next(r) 
        for row in r:
            if len(row) < 3:
                continue
            g_list.append(_parse_float(row[0]))
            gaia_list.append(_parse_float(row[1]))
            nancy_list.append(_parse_float(row[2]))
    return np.asarray(g_list), np.asarray(gaia_list), np.asarray(nancy_list)

def main() -> None:
    g, sigma_gaia, sigma_nancy = load_csv(CSV_PATH)
    order = np.argsort(g)
    g, sigma_gaia, sigma_nancy = g[order], sigma_gaia[order], sigma_nancy[order]

    finite = np.isfinite(sigma_gaia)
    g_fin, sig_fin = g[finite], sigma_gaia[finite]

    g_lo = float(g.min())
    g_hi = min(float(g.max()), G_MAX)
    in_win = g <= g_hi

    g_dense = np.linspace(g_lo, g_hi, N_CURVE)
    ok = g_dense <= g_fin[-1]
    curve_gaia = np.full_like(g_dense, np.inf)
    curve_gaia[ok] = np.interp(
        g_dense[ok], g_fin, sig_fin, left=sig_fin[0], right=sig_fin[-1]
    )
    curve_nancy = np.interp(
        g_dense, g, sigma_nancy, left=sigma_nancy[0], right=sigma_nancy[-1]
    )

    plug_ok = G_PLUG <= g_hi
    if plug_ok:
        pg = float(G_PLUG)
        if pg <= g_fin[-1]:
            v_gaia = float(
                np.interp(
                    pg, g_fin, sig_fin, left=sig_fin[0], right=sig_fin[-1]
                )
            )
        else:
            v_gaia = float("inf")
        v_nancy = float(
            np.interp(pg, g, sigma_nancy, left=sigma_nancy[0], right=sigma_nancy[-1])
        )
        print(f"G = {G_PLUG}:  Gaia DR4 alone = {v_gaia}  |  Gaia + NANCY = {v_nancy}")
    else:
        print(f"G_PLUG ({G_PLUG}) is above G_MAX ({G_MAX}); skipped.")

    fig, (ax0, ax1) = plt.subplots(2, 1, sharex=True, figsize=(7, 6), layout="constrained")
    show_gaia = np.isfinite(curve_gaia)
    ax0.plot(g_dense[show_gaia], curve_gaia[show_gaia], color="C0", lw=2, label="interp")
    ax0.scatter(
        g[finite & in_win],
        sigma_gaia[finite & in_win],
        c="k",
        s=25,
        zorder=5,
        label="table",
    )
    ax0.set_ylabel("Gaia DR4 alone")
    ax0.legend(loc="upper left", fontsize=9)
    ax0.grid(True, alpha=0.35)

    ax1.plot(g_dense, curve_nancy, color="C1", lw=2, label="interp")
    ax1.scatter(g[in_win], sigma_nancy[in_win], c="k", s=25, zorder=5, label="table")
    ax1.set_xlabel("Gaia G")
    ax1.set_ylabel("Gaia + NANCY")
    ax1.set_xlim(g_lo, g_hi)
    ax1.legend(loc="upper left", fontsize=9)
    ax1.grid(True, alpha=0.35)

    fig.suptitle(f"PM error vs Gaia G (G ≤ {g_hi:g})")
    plt.show()

if __name__ == "__main__":
    main()

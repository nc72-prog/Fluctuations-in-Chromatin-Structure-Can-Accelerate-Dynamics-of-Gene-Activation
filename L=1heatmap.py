import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

L = 1000
m = 500
u = 1e3
kon = 1e2

lambdas = [1e1, 1e2, 1e3]
koffs   = [u/lam**2 for lam in lambdas]

alphas = np.logspace(-2, 3, 50)
betas  = np.logspace(-2, 3, 50)

def compute_S1(L, m, u, koff):
    disc = np.sqrt((2*u + koff)**2 - 4*u**2)
    y = (2*u + koff - disc) / (2*u)
    if y <= 0 or y >= 1:
        return None
    logy = np.log(y)
    y2L    = np.exp(2*L*logy)
    y2m1   = np.exp((2*m-1)*logy)
    y2L12m = np.exp((2*L+1-2*m)*logy)
    S1 = (1+y)*(1-y2L) / ((1-y)*(1+y2m1)*(1+y2L12m))
    return S1

def compute_Tfree(L, m, u, kon, koff):
    S1 = compute_S1(L, m, u, koff)
    if S1 is None or S1 <= 0:
        return None
    return (1.0/kon)*(L/S1) + (1.0/koff)*((L-S1)/S1)

def compute_Tchr_L1(L, kon, alpha, beta):
    return L/kon + 1.0/alpha + beta/(alpha*kon) + (L-1)/beta


class CenterBandNorm(mcolors.Normalize):
    def __init__(self, vmin, vmax, lo=0.95, hi=1.05, clip=False):
        self.lo = lo
        self.hi = hi
        super().__init__(vmin, vmax, clip)

    def __call__(self, value, clip=None):
        value = np.ma.asarray(value, dtype=float)
        result = np.ma.zeros(value.shape)

        mask_low  = value <= self.lo
        mask_high = value >= self.hi
        mask_mid  = (~mask_low) & (~mask_high)

        result[mask_low]  = (value[mask_low] - self.vmin) / (self.lo - self.vmin) * 0.5
        result[mask_high] = 0.5 + (value[mask_high] - self.hi) / (self.vmax - self.hi) * 0.5
        result[mask_mid]  = 0.5

        return np.ma.clip(result, 0, 1)

sharp_rwg = mcolors.LinearSegmentedColormap.from_list(
    'sharp_rwg',
    [
        (0.00, 'darkred'),
        (0.25, 'red'),
        (0.50, 'white'),
        (0.75, 'green'),
        (1.00, 'darkgreen'),
    ]
)

n_panels = len(lambdas)
fig, axes = plt.subplots(1, n_panels, figsize=(4*n_panels, 4))
if n_panels == 1:
    axes = [axes]

labels = [rf'$\lambda=10^{{{int(np.log10(lam))}}}$' for lam in lambdas]

for idx, (lam, koff) in enumerate(zip(lambdas, koffs)):
    ratio = np.zeros((len(betas), len(alphas)))
    Tfree = compute_Tfree(L, m, u, kon, koff)

    for i, beta in enumerate(betas):
        for j, alpha in enumerate(alphas):
            Tchr = compute_Tchr_L1(L, kon, alpha, beta)
            if Tchr is not None and Tchr > 0:
                ratio[i, j] = Tfree / Tchr
            else:
                ratio[i, j] = np.nan

    vmax = max(np.nanmax(ratio), 1.06)
    vmin = min(np.nanmin(ratio), 0.94)
    norm = CenterBandNorm(vmin=vmin, vmax=vmax, lo=0.95, hi=1.05)

    ax = axes[idx]
    im = ax.pcolormesh(alphas, betas, ratio, norm=norm, cmap=sharp_rwg, shading='auto')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('α (s⁻¹)', fontsize=11)
    if idx == 0:
        ax.set_ylabel('β (s⁻¹)', fontsize=11)
    ax.set_title(f'{labels[idx]}, koff={koff:.0e}', fontsize=10)
    plt.colorbar(im, ax=ax, label='Tfree/Tchr')

plt.suptitle('Tfree / T_chr(L=1): green=chromatin faster, red=chromatin slower\nu=1e3, kon=1e2, L=1000, m=500', fontsize=12)
plt.tight_layout()
plt.savefig('heatmap_Lequals1.png', dpi=150, bbox_inches='tight')
plt.show()

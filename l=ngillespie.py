import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# --- fixed constants, must match the C++ file ---
L = 1000
m = 500
u = 1e3
koff = 10
kon = 1e2

def compute_S1(L, m, u, koff):
    disc = np.sqrt((2*u + koff)**2 - 4*u**2)
    y = (2*u + koff - disc) / (2*u)
    logy = np.log(y)
    y2L    = np.exp(2*L*logy)
    y2m1   = np.exp((2*m-1)*logy)
    y2L12m = np.exp((2*L+1-2*m)*logy)
    S1 = (1+y)*(1-y2L) / ((1-y)*(1+y2m1)*(1+y2L12m))
    return S1

def compute_Tfree(L, m, u, kon, koff):
    S1 = compute_S1(L, m, u, koff)
    return (1.0/kon)*(L/S1) + (1.0/koff)*((L-S1)/S1)

Tfree = compute_Tfree(L, m, u, kon, koff)

# --- load simulated T_chrom grid ---
data = np.loadtxt("gillespie_heatmap5.dat")
beta_col, alpha_col, T_sim_col = data[:, 0], data[:, 1], data[:, 2]

betas  = np.unique(beta_col)
alphas = np.unique(alpha_col)
nBeta, nAlpha = len(betas), len(alphas)

T_chrom = T_sim_col.reshape(nBeta, nAlpha)
ratio = Tfree / T_chrom


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
    [(0.00, 'darkred'), (0.25, 'red'), (0.50, 'white'), (0.75, 'green'), (1.00, 'darkgreen')]
)

vmax = max(np.nanmax(ratio), 1.06)
vmin = min(np.nanmin(ratio), 0.94)
norm = CenterBandNorm(vmin=vmin, vmax=vmax, lo=0.95, hi=1.05)

plt.figure(figsize=(7, 6))
im = plt.pcolormesh(alphas, betas, ratio, norm=norm, cmap=sharp_rwg, shading='auto')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\alpha$ (s$^{-1}$)', fontsize=20)
plt.ylabel(r'$\beta$ (s$^{-1}$)', fontsize=20)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.title(f'kon={kon}, koff={koff}', fontsize=18)
cbar = plt.colorbar(im)
cbar.set_label('Tfree / T00', fontsize=25)
plt.tight_layout()
plt.savefig('gillespie_heatmap.png', dpi=150)
plt.show()

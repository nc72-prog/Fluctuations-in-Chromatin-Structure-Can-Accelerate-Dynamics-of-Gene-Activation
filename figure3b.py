import numpy as np
import matplotlib.pyplot as plt

# ---- l = L model (whole-chain opening)
L, m = 1000, 500
kon, koff, u = 1e2, 1.0, 1e3             
alpha_values = [0.1, 1.0, 10.0]

def S0(beta):
    """Eq. 12 at s=0, written in an overflow-safe but algebraically identical form."""
    a = 2*u + koff + beta
    y = 2*u / (a + np.sqrt(a*a - 4*u*u))                    
    return (1+y)*(1 - y**(2*L)) / ((1-y)*(1 + y**(2*m-1))*(1 + y**(2*L+1-2*m)))

def T00_lL(beta, alpha):
    s = S0(beta)
    f = 1.0 + beta*L/alpha                                  
    T0 = f/kon*(L/s) + f/(beta + koff)*(L/s - 1.0)          
    return L/alpha + T0                                    

beta_ana = np.logspace(-3, 3, 300)

fig, ax = plt.subplots(figsize=(10, 7))
colors = ['tab:blue', 'tab:orange', 'tab:green']
for i, alpha in enumerate(alpha_values):
    ax.plot(beta_ana, T00_lL(beta_ana, alpha), color=colors[i], lw=4,
            label=rf"$\alpha$ = {alpha}")

ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"$\beta$ (s$^{-1}$)", fontsize=20)
ax.set_ylabel(r"$T_{00}^{(L)}$ (s)", fontsize=20)
ax.tick_params(labelsize=18); ax.legend(fontsize=18)
plt.tight_layout(); plt.savefig("gillespie_lL_beta.png", dpi=200); plt.show()

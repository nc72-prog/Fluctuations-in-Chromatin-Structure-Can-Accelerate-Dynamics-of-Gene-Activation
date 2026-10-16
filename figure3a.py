import numpy as np
import matplotlib.pyplot as plt
 
# ---- l = 1 model (single-site opening)
L   = 1000
kon = 1e2
alpha_values = [0.1, 1.0, 10.0]
 
def T00_l1(beta, alpha):
    tau = 1.0/(alpha*L) + 1.0/kon + beta/(kon*alpha*L)      
    return L*tau + (L - 1.0)/beta                           

beta_ana = np.logspace(-2, 5, 300)
 
fig, ax = plt.subplots(figsize=(10, 7))
colors = ['tab:blue', 'tab:orange', 'tab:green']
for i, alpha in enumerate(alpha_values):
    ax.plot(beta_ana, T00_l1(beta_ana, alpha), color=colors[i], lw=4,
            label=rf"$\alpha$ = {alpha}")
 
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"$\beta$ (s$^{-1}$)", fontsize=20)
ax.set_ylabel(r"$T_{00}^{(1)}$ (s)", fontsize=20)
ax.tick_params(labelsize=18); ax.legend(fontsize=18)
plt.tight_layout(); plt.savefig("analytical_l1_beta.png", dpi=200); plt.show()
 

import numpy as np
import matplotlib.pyplot as plt

def plot_L2_convergence(h, L2, output_path="plotL2.png"):
    plt.figure(figsize=(6, 4))
    plt.loglog(h, L2, 'bo-', label='L2 error')
    plt.xlabel('Element size (h)')
    plt.ylabel('L2 error norm')
    plt.title('L2 Error Norm vs Element Size')
    plt.grid(True, which="both", ls="--")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.show()

    log_h = np.log(h)
    log_L2 = np.log(L2)
    slope, _ = np.polyfit(log_h, log_L2, 1)
    print(f"Convergence rate (slope) = {slope:.4f}")
    return slope

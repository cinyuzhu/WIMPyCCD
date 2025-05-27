import numpy as np
import matplotlib.pyplot as plt

# Define conversion factor
def plot_differential_rate(E_Rs, dR_values, labels, 
                           ref_index = 1, 
                           panel_type='ratio',  # 'residual', 'weighted_residual', or 'ratio'
                           show = False):
    """
    Plots differential rate curves with an additional comparison panel.
    """
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), gridspec_kw={'height_ratios': [3, 1]}, sharex=True)

    dR_values = np.array(dR_values)
    # Main plot
    for i, (dR, label) in enumerate(zip(dR_values, labels)):
            if i == ref_index:
                axes[0].plot(E_Rs, dR, linestyle=(0,(5,1)), linewidth=2.5, label=f"{label} (Reference)")
            else:
                axes[0].plot(E_Rs, dR, label=label)
    axes[0].set_ylabel("dR/dE (day⁻¹ kg⁻¹ KeV⁻¹)")
    axes[0].set_title("Differential Rate")
    axes[0].legend()

    # Comparison panel 
    dR_reference = dR_values[ref_index]
    if panel_type == 'ratio':
        comp_values = [dR / dR_reference for dR in dR_values]
        ylabel = "Ratio"
        ref_line = 1
    elif panel_type == 'residual':
        comp_values = [dR - dR_reference for dR in dR_values]
        ylabel = "Residual"
        ref_line = 0
    elif panel_type == 'weighted_residual':
        # Avoid division by zero — assume Poisson for now
        comp_values = [(dR - dR_reference) / np.sqrt(dR_reference + 1e-12) for dR in dR_values]
        ylabel = "weighted_residual"
        ref_line = 0
    else:
        raise ValueError(f"Unknown comparison_type: {panel_type}")
    # Compute ratio relative to reference
    label_ref = labels[ref_index]
    for i, (comp, label) in enumerate(zip(comp_values, labels)):
        if i == ref_index:
            axes[1].plot(
                E_Rs, comp, linestyle=(0, (5, 1)), linewidth=2.5,
                label=f"{label} (Reference)"
            )
        else:
            if panel_type == "ratio":
                # comp_label = f"\[{label}\] / \[{label_ref}\]"
                comp_label = f"{label} / {label_ref}"
            elif panel_type == "residual":
                comp_label = f"{label} - {label_ref}"
            elif panel_type == "weighted_residual":
                comp_label = f"({label} - {label_ref}) / √{label_ref}"
            else:
                comp_label = label  # fallback

            axes[1].plot(E_Rs, comp, label=comp_label)

    axes[1].axhline(ref_line, color="gray", linestyle="--")
    axes[1].set_xlabel("Energy (KeV)")
    axes[1].set_ylabel(ylabel)
    axes[1].legend()

    plt.tight_layout()
    if show:
        plt.show()

    return fig, axes

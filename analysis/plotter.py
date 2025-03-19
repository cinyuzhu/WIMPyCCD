import numpy as np
import matplotlib.pyplot as plt

# Define conversion factor


def plot_differential_rate_with_ratio(E_Rs_GeV, dR_values, labels, 
                                      unit_conversion=True, reference_index = 1,
                                      xlabel="Energy (KeV)", ylabel="dR/dE (day⁻¹ kg⁻¹ KeV⁻¹)", 
                                      title="Differential Rate", ratio_rage= [0.5, 1.5]):
    """
    Plots differential rate curves with an additional ratio panel.

    Parameters:
    - E_Rs_GeV: Array of energy values in GeV
    - dR_values: arrays representing differential rates
    - labels: List of strings, corresponding to the curves
    - reference_index: Index of the reference curve for ratio calculation (the second curve by default)
    - xlabel, ylabel, title: Plot labels for customization
    """
    if unit_conversion:
        GeV_to_KeV = 1e6  # 1 GeV = 10^6 KeV
    else:
        GeV_to_KeV = 1

    # Convert energy to KeV
    E_Rs_KeV = E_Rs_GeV * GeV_to_KeV

    # Convert differential rates to KeV⁻¹
    dR_values_KeV = [np.array(dR) / GeV_to_KeV for dR in dR_values]

    # Choose reference curve for ratio calculation (default: second curve, dRs1)
    dR_reference = dR_values_KeV[reference_index]

    # Compute ratio relative to reference
    ratio_values = [dR / dR_reference for dR in dR_values_KeV]

    # Create figure with two subplots (main plot + ratio panel)
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), gridspec_kw={'height_ratios': [3, 1]}, sharex=True)

    # Main plot
    for i, (dR, label) in enumerate(zip(dR_values_KeV, labels)):
            if i == reference_index:
                axes[0].plot(E_Rs_KeV, dR, linestyle=(0,(5,1)), linewidth=4, label=f"{label} (Reference)")
            else:
                axes[0].plot(E_Rs_KeV, dR, label=label)
    axes[0].set_ylabel(ylabel)
    axes[0].set_title(title)
    axes[0].legend()

    # Ratio panel
    for i, (ratio, label) in enumerate(zip(ratio_values, labels)):
        if i == reference_index:
            axes[1].plot(E_Rs_KeV, ratio, linestyle=(0,(5,1)), linewidth=4, label=f"{label} (Reference)")
        else:
            axes[1].plot(E_Rs_KeV, ratio, label=label)

    axes[1].axhline(1, color="gray", linestyle="--")  # Reference line at y=1
    axes[1].set_xlabel(xlabel)
    axes[1].set_ylabel("Ratio")
    axes[1].set_ylim(ratio_rage)


    plt.tight_layout()
    plt.show()

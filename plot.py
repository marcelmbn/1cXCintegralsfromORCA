"""
Plot the 1c XC integral data.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def plot_onexc_ints(array: np.ndarray) -> None:
    """
    Plot the 1c XC integral data as a multi-line chart using seaborn.

    The shape corresponds to this initialization: onecxcints = np.zeros((9, 104))
    - array.shape[0] corresponds to different shell interactions (1..9).
    - array.shape[1] corresponds to different elements (1..104).

    Args:
        array (np.ndarray): The 1c XC integral data with shape (9, 104).

    Returns:
        None
    """
    # Optional: set a nice theme
    sns.set_theme(style="whitegrid", context="talk")

    # Build a DataFrame from the array for easier plotting with seaborn
    # We’ll label each row (shell_interaction) and each column (element)
    # so that we can use seaborn's `lineplot` with a "hue" grouping.
    data = []
    for shell_idx in range(array.shape[0]):
        for elem_idx in range(array.shape[1]):
            shell_interaction: str
            if shell_idx == 0:
                shell_interaction = "s - p"
            elif shell_idx == 1:
                shell_interaction = "p - p'"
            elif shell_idx == 2:
                shell_interaction = "s - d"
            elif shell_idx == 3:
                shell_interaction = "p - d"
            elif shell_idx == 4:
                shell_interaction = "d - d'"
            elif shell_idx == 5:
                shell_interaction = "s - f"
            elif shell_idx == 6:
                shell_interaction = "p - f"
            elif shell_idx == 7:
                shell_interaction = "d - f"
            elif shell_idx == 8:
                shell_interaction = "f - f'"
            data.append(
                {
                    "Shell Interaction": shell_interaction,
                    "Element Number": elem_idx + 1,
                    "Integral Value": array[shell_idx, elem_idx],
                }
            )

    df = pd.DataFrame(data)

    plt.figure(figsize=(10, 6))
    sns.lineplot(
        data=df,
        x="Element Number",
        y="Integral Value",
        hue="Shell Interaction",
        palette="tab10",  # Choose a nice color palette
    )

    plt.title("shell-averaged one-center exchange integrals")
    plt.xlabel("element number")
    plt.ylabel("integral value / a.u.")
    plt.legend(title="shell interaction", loc="best")
    plt.tight_layout()
    plt.savefig("onecxcints.png", dpi=300)
    plt.show()

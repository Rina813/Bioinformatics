"""
Command-line interface for qPCR Analyzer.

Run:
    python scripts/qpcr_cli.py
"""

import matplotlib.pyplot as plt

from my_module.functions import summarize_target_from_file


def plot_cq_summary(summary_df, target_name):
    """
    Plot a bar graph of mean Cq per base_sample with error bars (std).
    """
    if summary_df.empty:
        print("Nothing to plot (summary is empty).")
        return

    summary_df = summary_df.sort_values("base_sample")

    x_labels = summary_df["base_sample"].tolist()
    mean_cq = summary_df["mean_cq"].tolist()
    std_cq = summary_df["std_cq"].tolist()

    x_positions = range(len(x_labels))

    plt.figure(figsize=(9, 5))
    plt.bar(x_positions, mean_cq, yerr=std_cq, capsize=5)
    plt.xticks(x_positions, x_labels, rotation=45, ha="right")
    plt.ylabel("Cq")
    plt.xlabel("Sample")
    plt.title(f"Mean Cq per Sample for Target '{target_name}'")
    plt.tight_layout()
    plt.show()


def main():
    print("=" * 60)
    print(" qPCR Analyzer: Mean Cq per Sample (duplicates merged) ")
    print("=" * 60)

    filename = input("Enter CSV filename (e.g., 20251103_1.csv):\n> ").strip()
    target = input("\nEnter Target (gene) name (e.g., Actb):\n> ").strip()

    try:
        summary = summarize_target_from_file(filename, target)
    except FileNotFoundError:
        print(f"\n[Error] File not found: {filename}")
        print("Make sure the CSV is uploaded into this folder, or provide a relative path.")
        return
    except KeyError as e:
        print(f"\n[Error] {e}")
        print("Your CSV must contain columns: Target, Sample, Cq")
        return
    except Exception as e:
        print("\n[Error] Unexpected failure.")
        print(f"Details: {e}")
        return

    if summary.empty:
        print(f"\nNo rows found for Target = '{target}'. Check spelling/case.")
        return

    print("\nSummary:")
    print(summary.to_string(index=False))

    do_plot = input("\nShow plot? (y/n): ").strip().lower()
    if do_plot == "y":
        plot_cq_summary(summary, target)

    do_save = input("\nSave summary to CSV? (y/n): ").strip().lower()
    if do_save == "y":
        out_name = input("Output filename (e.g., Actb_summary.csv):\n> ").strip()
        if not out_name:
            out_name = f"{target}_summary.csv"
        summary.to_csv(out_name, index=False)
        print(f"Saved: {out_name}")

    print("\nDone!")


if __name__ == "__main__":
    main()

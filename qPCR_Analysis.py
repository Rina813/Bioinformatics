import csv
import math
filename = r"/Users/rinamatsuo/Documents/Karin_lab/2025_11/20251105/20251103_1.csv"


def read_qpcr_csv(filename):
    """Read qPCR CSV into a list of dicts."""
    rows = []
    with open(filename, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            row["ct"] = float(row["ct"])
            rows.append(row)
    return rows


def group_ct_values(data):
    """
    Group Ct values by (condition, gene, sample_id).

    Returns:
        dict[(condition, gene)][sample_id] = ct
    """
    grouped = {}
    for row in data:
        cond = row["condition"]
        gene = row["gene"]
        sample = row["sample_id"]
        ct = row["ct"]

        key = (cond, gene)
        if key not in grouped:
            grouped[key] = {}
        grouped[key][sample] = ct
    return grouped


def compute_delta_ct(grouped, housekeeping_gene):
    """
    Compute ΔCt for each (condition, gene, sample).

    Returns:
        delta_ct[(condition, gene)] = list of ΔCt values
    """
    delta_ct = {}

    # First get housekeeping Ct per (condition, sample)
    housekeeping = {}
    for (cond, gene), samples in grouped.items():
        if gene == housekeeping_gene:
            for sample_id, ct in samples.items():
                housekeeping[(cond, sample_id)] = ct

    # Now compute ΔCt for each gene using housekeeping
    for (cond, gene), samples in grouped.items():
        if gene == housekeeping_gene:
            # skip housekeeping itself (ΔCt for housekeeping would be 0)
            continue

        key = (cond, gene)
        delta_ct[key] = []

        for sample_id, ct in samples.items():
            hk_key = (cond, sample_id)
            if hk_key not in housekeeping:
                # housekeeping Ct missing for this sample/condition
                # you could print a warning or skip
                continue

            hk_ct = housekeeping[hk_key]
            dct = ct - hk_ct
            delta_ct[key].append(dct)

    return delta_ct


def average_delta_ct(delta_ct_dict):
    """Average ΔCt values per (condition, gene)."""
    mean_delta_ct = {}
    for key, values in delta_ct_dict.items():
        if len(values) == 0:
            continue
        mean_delta_ct[key] = sum(values) / len(values)
    return mean_delta_ct


def compute_ddct(mean_delta_ct, control_condition):
    """
    Compute ΔΔCt compared to control_condition.

    Returns:
        ddct[(condition, gene)] = ΔΔCt
    """
    ddct = {}

    # Find all genes
    genes = set(g for (cond, g) in mean_delta_ct.keys())

    for gene in genes:
        control_key = (control_condition, gene)
        if control_key not in mean_delta_ct:
            # no control data for this gene
            continue

        control_dct = mean_delta_ct[control_key]

        for (cond, g) in mean_delta_ct.keys():
            if g != gene:
                continue
            dct = mean_delta_ct[(cond, g)]
            ddct[(cond, g)] = dct - control_dct

    return ddct


def compute_fold_change(ddct_dict):
    """Compute fold change = 2^(-ΔΔCt) for each (condition, gene)."""
    fold_change = {}
    for key, value in ddct_dict.items():
        fold_change[key] = 2 ** (-value)
    return fold_change


def write_results_to_csv(filename, mean_delta_ct, ddct, fold_change):
    """Write summary results to a CSV file."""
    with open(filename, "w", newline="") as f:
        fieldnames = ["gene", "condition", "mean_delta_ct", "ddct", "fold_change"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        # sort keys for nice output
        all_keys = sorted(mean_delta_ct.keys())

        for (cond, gene) in all_keys:
            row = {
                "gene": gene,
                "condition": cond,
                "mean_delta_ct": round(mean_delta_ct.get((cond, gene), float("nan")), 4),
                "ddct": round(ddct.get((cond, gene), float("nan")), 4)
                if (cond, gene) in ddct else "",
                "fold_change": round(fold_change.get((cond, gene), float("nan")), 4)
                if (cond, gene) in fold_change else "",
            }
            writer.writerow(row)


def main():
    print("Welcome to qPCR Analyzer (ΔΔCt calculator)")
    input_file = input("Enter input CSV filename (e.g., data.csv): ")
    hk_gene = input("Enter housekeeping gene name (e.g., GAPDH): ")
    control_cond = input("Enter control condition name (e.g., control): ")
    output_file = input("Enter output CSV filename (e.g., results.csv): ")

    data = read_qpcr_csv(input_file)
    grouped = group_ct_values(data)
    delta_ct = compute_delta_ct(grouped, hk_gene)
    mean_delta_ct = average_delta_ct(delta_ct)
    ddct = compute_ddct(mean_delta_ct, control_cond)
    fold_change = compute_fold_change(ddct)

    write_results_to_csv(output_file, mean_delta_ct, ddct, fold_change)

    print("\nAnalysis complete!")
    print(f"Results saved to {output_file}")
    print("Fold changes are relative to the control condition using 2^(-ΔΔCt).")


if __name__ == "__main__":
    main()

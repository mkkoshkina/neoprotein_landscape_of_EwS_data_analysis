"""
Patient simulation pipeline for HLA binder counts from ORF FASTA.
 
Steps:
1. Parse input FASTA
2. Sliding window 8-11 aa → write .pep files → run NetMHCpan
3. Parse NetMHCpan outputs → build per-allele SB/WB/B dicts
4. Simulate N patients by sampling HLA alleles from population frequencies
5. Count binders per patient and plot histogram
"""
 
import os
import random
import re
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
 
 
# ─────────────────────────────────────────────
# USER SETTINGS
# ─────────────────────────────────────────────
 
INPUT_FASTA    = "orfs.fasta"
NETMHCPAN_PATH = "/home/mariia/Desktop/netMHCpan-4.1/netMHCpan"
OUTPUT_DIR     = "netmhcpan_outputs_simulation"
RESULTS_CSV    = "simulation_results.csv"
PLOT_OUT       = "simulation_histogram.pdf"
 
N_PATIENTS   = 10000
RANDOM_SEED  = 100
SB_THRESHOLD = 0.5
WB_THRESHOLD = 2.0
MHC_LENGTHS  = [8, 9, 10, 11]
 
NETMHCPAN_ALLELES = [
    "HLA-A01:01", "HLA-A02:01", "HLA-A03:01", "HLA-A11:01",
    "HLA-A24:02", "HLA-A29:02",
    "HLA-B07:02", "HLA-B08:01", "HLA-B15:01", "HLA-B18:01",
    "HLA-B44:03", "HLA-B51:01",
    "HLA-C02:02", "HLA-C03:03", "HLA-C04:01", "HLA-C05:01",
    "HLA-C07:01", "HLA-C07:02",
]
 
A_ALLELES_PS = {
    "A0101": 0.1301, "A0201": 0.2901, "A0301": 0.1226,
    "A1101": 0.05,   "A2402": 0.0935, "A2902": 0.0566,
    "other": 0.2571,
}
B_ALLELES_PS = {
    "B0702": 0.1137, "B0801": 0.0953, "B1501": 0.0536,
    "B1801": 0.058,  "B4403": 0.0801, "B5101": 0.0774,
    "other": 0.5755,
}
C_ALLELES_PS = {
    "C0202": 0.0672, "C0303": 0.0566, "C0401": 0.1339,
    "C0501": 0.1573, "C0701": 0.1814, "C0702": 0.0944,
    "other": 0.3092,
}
 
 
# ─────────────────────────────────────────────
# UTILITIES
# ─────────────────────────────────────────────
 
def safe_id(prot_id: str) -> str:
    return re.sub(r"[|/\\:*?\"<> ]", "_", prot_id)
 
def hla_long_to_short(string: str) -> str:
    return string.replace("HLA-", "").replace(":", "")
 
 
# ─────────────────────────────────────────────
# STEP 1: Load FASTA
# ─────────────────────────────────────────────
 
def load_fasta(path: str) -> dict:
    records = {record.id: str(record.seq) for record in SeqIO.parse(path, "fasta")}
    print(f"Loaded {len(records)} ORFs from {path}")
    return records
 
 
# ─────────────────────────────────────────────
# STEP 2: Sliding window → .pep files → run NetMHCpan
# ─────────────────────────────────────────────
 
def write_pep_file(sequence: str, path: str) -> None:
    with open(path, "w") as f:
        for l in MHC_LENGTHS:
            for i in range(len(sequence) - l + 1):
                f.write(sequence[i:i+l] + "\n")
 
 
def run_netmhcpan(all_proteins: dict) -> None:
    alleles_str = ",".join(NETMHCPAN_ALLELES)
    commands = []
    for prot_id, sequence in all_proteins.items():
        sid = safe_id(prot_id)
        pep_file = os.path.join(OUTPUT_DIR, f"{sid}.pep")
        xls_file = os.path.join(OUTPUT_DIR, f"{sid}_NetMHCpan_out.xls")
        write_pep_file(sequence, pep_file)
        commands.append(f"{NETMHCPAN_PATH} -a {alleles_str} -p {pep_file} -xls -xlsfile {xls_file}\n")
 
    sh_file = "netmhcpan_commands_simulation.sh"
    with open(sh_file, "w") as f:
        f.write("#!/bin/bash\n\n")
        f.writelines(commands)
 
    print(f"Running NetMHCpan on {len(all_proteins)} ORFs...")
    subprocess.run(["bash", sh_file], check=True)
    print("NetMHCpan done.")
 
 
# ─────────────────────────────────────────────
# STEP 3: Parse outputs → per-allele dicts
# ─────────────────────────────────────────────
 
def parse_netmhcpan_output(xls_file: str) -> dict:
    with open(xls_file) as f:
        hlas_long = [h.strip() for h in f.readline().rstrip("\n").split("\t\t\t\t") if h.strip()]
 
    df = pd.read_csv(xls_file, index_col=1, sep="\t", header=1)
    df = df.rename(columns={"EL_Rank": "EL_Rank.0"}).reset_index()
    df = df[df["NB"] != 0]
 
    el_rank_cols = [c for c in df.columns if c.startswith("EL_Rank")]
    df = df.rename(columns={el_rank_cols[i]: hlas_long[i] for i in range(min(len(el_rank_cols), len(hlas_long)))})
 
    result = {}
    for hla_long in hlas_long:
        if hla_long not in df.columns:
            continue
        hla_short = hla_long_to_short(hla_long)
        sb = len(df[df[hla_long] < SB_THRESHOLD])
        wb = len(df[(df[hla_long] >= SB_THRESHOLD) & (df[hla_long] < WB_THRESHOLD)])
        result[hla_short] = (sb, wb, sb + wb)
    return result
 
 
def build_allele_dicts(all_proteins: dict) -> tuple:
    SB = {hla_long_to_short(h): 0 for h in NETMHCPAN_ALLELES}
    WB = {hla_long_to_short(h): 0 for h in NETMHCPAN_ALLELES}
    B  = {hla_long_to_short(h): 0 for h in NETMHCPAN_ALLELES}
 
    for prot_id in all_proteins:
        xls_file = os.path.join(OUTPUT_DIR, f"{safe_id(prot_id)}_NetMHCpan_out.xls")
        try:
            counts = parse_netmhcpan_output(xls_file)
            for hla_short, (sb, wb, b) in counts.items():
                SB[hla_short] += sb
                WB[hla_short] += wb
                B[hla_short]  += b
        except Exception as e:
            print(f"Warning: could not parse {xls_file}: {e}")
 
    for d in (SB, WB, B):
        d["other"] = 0
 
    return SB, WB, B
 
 
# ─────────────────────────────────────────────
# STEP 4: Patient simulation
# ─────────────────────────────────────────────
 
def sample_patients(n: int) -> pd.DataFrame:
    random.seed(RANDOM_SEED)
    rows = []
    for i in range(n):
        hlas = [
            random.choices(list(A_ALLELES_PS.keys()), weights=list(A_ALLELES_PS.values()))[0],
            random.choices(list(A_ALLELES_PS.keys()), weights=list(A_ALLELES_PS.values()))[0],
            random.choices(list(B_ALLELES_PS.keys()), weights=list(B_ALLELES_PS.values()))[0],
            random.choices(list(B_ALLELES_PS.keys()), weights=list(B_ALLELES_PS.values()))[0],
            random.choices(list(C_ALLELES_PS.keys()), weights=list(C_ALLELES_PS.values()))[0],
            random.choices(list(C_ALLELES_PS.keys()), weights=list(C_ALLELES_PS.values()))[0],
        ]
        rows.append({"patient_n": f"patient_{i+1}", "HLAs": hlas})
    return pd.DataFrame(rows)
 
 
def count_binders(hlas: list, SB: dict, WB: dict, B: dict) -> tuple:
    return (
        sum(SB.get(h, 0) for h in hlas),
        sum(WB.get(h, 0) for h in hlas),
        sum(B.get(h, 0)  for h in hlas),
    )
 
 
# ─────────────────────────────────────────────
# STEP 5: Plot
# ─────────────────────────────────────────────
 
def plot_histogram(df: pd.DataFrame) -> None:
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.2)
 
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.histplot(data=df, x="SB", stat="probability", color="#1b365d", ax=ax)
    ax.set_xlabel("Strong binder count per patient", fontsize=12)
    ax.set_ylabel("Probability", fontsize=12)
    ax.set_title(
        f"Distribution of strong binder counts in {N_PATIENTS:,} simulated patients\n"
        f"mean = {df['SB'].mean():.2f}, std = {df['SB'].std():.2f}",
        fontsize=12,
    )
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    plt.savefig(PLOT_OUT, bbox_inches="tight", dpi=300)
    print(f"Plot saved to {PLOT_OUT}")
    plt.show()
 
 
# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
 
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
 
    all_proteins = load_fasta(INPUT_FASTA)
    run_netmhcpan(all_proteins)
 
    SB, WB, B = build_allele_dicts(all_proteins)
 
    df = sample_patients(N_PATIENTS)
    df[["SB", "WB", "B"]] = pd.DataFrame(
        df["HLAs"].apply(lambda h: count_binders(h, SB, WB, B)).tolist(),
        index=df.index,
    )
 
    df.to_csv(RESULTS_CSV, index=False)
    print(f"Results saved to {RESULTS_CSV}")
    plot_histogram(df)
 
 
if __name__ == "__main__":
    main()
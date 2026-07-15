"""
Pipeline to calculate tryptic peptide density and HLA binder density
for neogenes vs canonical proteins. v2 — fast regex trypsin, main() guard.

Steps:
1. Parse two FASTAs (neogenes + canonical)
2. Tryptic digest via regex → filter 7-25 aa → compute density
3. Sliding window 8-11 aa → write .pep files → print NetMHCpan commands
4. Parse NetMHCpan outputs → count SB → compute density
5. Output a dataframe ready to plot

Run once to generate .pep files and print NetMHCpan commands.
After running NetMHCpan, run again — step 4 will pick up the outputs.
"""

import os
import re
import pandas as pd
import subprocess
from Bio import SeqIO


# ─────────────────────────────────────────────
# USER SETTINGS
# ─────────────────────────────────────────────

NEOGENE_FASTA  = "neogenes.fasta"
CANONICAL_FASTA = "canonical.fasta"
N_PROTEINS = 364

NETMHCPAN_PATH = "/home/mariia/Desktop/netMHCpan-4.1/netMHCpan"
HLA_ALLELES = "HLA-B35:03,HLA-A30:04,HLA-B15:03,HLA-B44:02,HLA-C04:01,HLA-C02:10,HLA-A03:01,HLA-C12:02,HLA-B07:02,HLA-B51:01,HLA-A02:01,HLA-A30:01,HLA-B40:02,HLA-C16:01,HLA-C16:02,HLA-B35:01,HLA-A01:01,HLA-C05:01,HLA-C06:02,HLA-C02:02,HLA-A33:03,HLA-B15:01,HLA-C15:02,HLA-B51:08,HLA-B57:01,HLA-B52:01,HLA-C03:04,HLA-C07:02,HLA-A32:01,HLA-B44:03,HLA-A68:01,HLA-C03:03,HLA-B50:01"  # replace with cohort alleles

OUTPUT_DIR  = "netmhcpan_outputs"
RESULTS_CSV = "detectability_results.csv"

MISSED_CLEAVAGES = 2
MIN_PEP_LEN = 7
MAX_PEP_LEN = 25
MHC_LENGTHS = [8, 9, 10, 11]
SB_THRESHOLD = 0.5  # EL_Rank cutoff for strong binders


# ─────────────────────────────────────────────
# STEP 1: Load FASTAs
# ─────────────────────────────────────────────

def load_fasta(path: str, n: int, label: str) -> dict:
    records = {}
    for record in SeqIO.parse(path, "fasta"):
        records[record.id] = str(record.seq)
        if len(records) >= n:
            break
    print(f"Loaded {len(records)} {label} proteins from {path}")
    return records


# ─────────────────────────────────────────────
# STEP 2: Tryptic digest → density (fast regex)
# ─────────────────────────────────────────────

# Trypsin: cleave after K or R, not before P
_TRYPSIN_RE = re.compile(r"(?<=[KR])(?!P)")

def tryptic_density(
    sequence: str,
    missed_cleavages: int = MISSED_CLEAVAGES,
    min_len: int = MIN_PEP_LEN,
    max_len: int = MAX_PEP_LEN,
) -> float:
    """
    Fast tryptic density via regex cleavage + index arithmetic.
    No string allocation for out-of-range peptides.
    """
    if not sequence:
        return 0.0

    sites = [0] + [m.end() for m in _TRYPSIN_RE.finditer(sequence)] + [len(sequence)]
    n = len(sites) - 1  # number of base fragments

    count = 0
    for i in range(n):
        for j in range(i + 1, min(i + missed_cleavages + 2, n + 1)):
            pep_len = sites[j] - sites[i]
            if pep_len > max_len:
                break  # fragments only get longer from here
            if pep_len >= min_len:
                count += 1

    return count / len(sequence)


# ─────────────────────────────────────────────
# STEP 3: Sliding window → .pep files + NetMHCpan commands
# ─────────────────────────────────────────────

def write_pep_file(sequence: str, path: str, lengths: list = MHC_LENGTHS) -> None:
    with open(path, "w") as f:
        for l in lengths:
            for i in range(len(sequence) - l + 1):
                f.write(sequence[i:i+l] + "\n")


def write_netmhcpan_commands(all_proteins: dict, out_file: str = "netmhcpan_commands.sh") -> None:
    with open(out_file, "w") as f:
        f.write("#!/bin/bash\n\n")
        for prot_id in all_proteins:
            pep_file = os.path.join(OUTPUT_DIR, f"{prot_id}.pep")
            xls_file = os.path.join(OUTPUT_DIR, f"{prot_id}_NetMHCpan_out.xls")
            f.write(f"{NETMHCPAN_PATH} -a {HLA_ALLELES} -p {pep_file} -xls -xlsfile {xls_file}\n")
    print(f"NetMHCpan commands written to {out_file}")
    print(f"Run with: bash {out_file}")


# ─────────────────────────────────────────────
# STEP 4: Parse NetMHCpan outputs → SB density
# ─────────────────────────────────────────────

def parse_netmhcpan_output(xls_file: str) -> int:
    """
    Returns total SB count across all alleles.
    Same peptide binding N alleles counts N times.
    """
    with open(xls_file) as f:
        hlas = [h.strip() for h in f.readline().rstrip("\n").split("\t\t\t\t") if h.strip()]

    df = pd.read_csv(xls_file, index_col=1, sep="\t", header=1)
    df = df.rename(columns={"EL_Rank": "EL_Rank.0"}).reset_index()
    df = df[df["NB"] != 0]

    el_rank_cols = [c for c in df.columns if c.startswith("EL_Rank")]
    df = df.rename(columns={el_rank_cols[i]: hlas[i] for i in range(min(len(el_rank_cols), len(hlas)))})

    return sum(
        len(df[df[hla] < SB_THRESHOLD])
        for hla in hlas
        if hla in df.columns
    )


# ─────────────────────────────────────────────
# STEP 5: Assemble results dataframe
# ─────────────────────────────────────────────

def build_results(all_proteins: dict) -> pd.DataFrame:
    rows = []
    for prot_id, sequence in all_proteins.items():
        label = "neogene" if prot_id.startswith("NEO_") else "canonical"
        xls_file = os.path.join(OUTPUT_DIR, f"{prot_id}_NetMHCpan_out.xls")

        tryp_dens = tryptic_density(sequence)

        sb_count, hla_dens = None, None
        if os.path.exists(xls_file):
            try:
                sb_count = parse_netmhcpan_output(xls_file)
                hla_dens = sb_count / len(sequence) if sequence else 0.0
            except Exception as e:
                print(f"Warning: could not parse {xls_file}: {e}")

        rows.append({
            "protein_id":              prot_id,
            "label":                   label,
            "length":                  len(sequence),
            "tryptic_peptide_density": tryp_dens,
            "hla_sb_count":            sb_count,
            "hla_sb_density":          hla_dens,
        })

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    neogenes   = load_fasta(NEOGENE_FASTA,   N_PROTEINS, "neogene")
    canonicals = load_fasta(CANONICAL_FASTA, N_PROTEINS, "canonical")

    all_proteins = {
        **{f"NEO_{k}": v for k, v in neogenes.items()},
        **{f"CAN_{k}": v for k, v in canonicals.items()},
    }

    # Write .pep files and print commands
    for prot_id, sequence in all_proteins.items():
        write_pep_file(sequence, os.path.join(OUTPUT_DIR, f"{prot_id}.pep"))
    write_netmhcpan_commands(all_proteins)
    subprocess.run(["bash", "netmhcpan_commands.sh"], check=True)

    # Parse whatever NetMHCpan outputs exist and save results
    df = build_results(all_proteins)
    df.to_csv(RESULTS_CSV, index=False)
    print(f"\nResults saved to {RESULTS_CSV}")

    n_parsed = df["hla_sb_density"].notna().sum()
    print(f"{n_parsed}/{len(df)} proteins have NetMHCpan output parsed.")
    if n_parsed > 0:
        print(df.groupby("label")[["tryptic_peptide_density", "hla_sb_density"]].describe())


if __name__ == "__main__":
    main()

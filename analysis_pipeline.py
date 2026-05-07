#!/usr/bin/env python3


import os
import sys
import argparse
import subprocess
from shutil import which
import math
import textwrap

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib import distances


# ----------------------------- Utility functions -----------------------------

def shutil_which(name):
    return which(name)


def run_freesasa(pdb_path, out_txt):
    """Run freesasa on pdb_path and write output to out_txt. Return True on success."""
    cmd = ["freesasa", pdb_path, "-o", out_txt]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return True
    except Exception as e:
        return False


def parse_freesasa_total(txt_path):
    try:
        with open(txt_path) as fh:
            for line in fh:
                if line.strip().upper().startswith("TOTAL"):
                    parts = line.split()
                    return float(parts[-1])
    except Exception:
        return None
    return None


# ----------------------------- Analysis helpers -----------------------------

#def residue_contacts(target_atoms, binder_atoms, cutoff=4.0, add_target_res_offset=0):
   # """Return sorted list of contacting residue pairs (target_resid, target_resname), (binder_resid, binder_resname).
   # target_atoms and binder_atoms are MDAnalysis AtomGroup (non-hydrogen recommended).
 #   add_target_res_offset: if your PDB target residue numbering was shifted and you want to add an offset (e.g. +23).
  #  """
  #  if len(target_atoms) == 0 or len(binder_atoms) == 0:
    #    return []
   # d = distances.distance_array(target_atoms.positions, binder_atoms.positions)
    #i, j = np.where(d <= cutoff)
    #pairs = set()
    #for a, b in zip(i, j):
        #resA = (int(target_atoms[a].resid + add_target_res_offset), target_atoms[a].resname)
       # resB = (int(binder_atoms[b].resid), binder_atoms[b].resname)
       # pairs.add((resA, resB))
    #return sorted(pairs, key=lambda x: (x[0][0], x[1][0]))

def residue_contacts(target_atoms, binder_atoms, cutoff=4.0, add_target_res_offset=0):
    """
    Returns:
    - contacts: sorted list of ((resA_id, resA_name), (resB_id, resB_name))
    - hydrophobic_patches: list of contacts specifically between two hydrophobic residues
    """
    # 1. Define Hydrophobic Lookup Table
    HYDROPHOBIC_AAS = {'VAL', 'ILE', 'LEU', 'PHE', 'MET', 'TRP', 'ALA', 'PRO'}

    if len(target_atoms) == 0 or len(binder_atoms) == 0:
        return [], []

    # 2. Calculate Distances
    d = distances.distance_array(target_atoms.positions, binder_atoms.positions)
    i, j = np.where(d <= cutoff)
    
    pairs = set()
    hydrophobic_patches = []

    # 3. Analyze Contacts
    for a, b in zip(i, j):
        # Extract residue info
        resA_id = int(target_atoms[a].resid + add_target_res_offset)
        resA_name = target_atoms[a].resname
        
        resB_id = int(binder_atoms[b].resid)
        resB_name = binder_atoms[b].resname
        
        resA_info = (resA_id, resA_name)
        resB_info = (resB_id, resB_name)

        
        # Add to unique set of pairs
        if (resA_info, resB_info) not in pairs:
            pairs.add((resA_info, resB_info))
            #pairs.add(dd)
            
            # 4. Identify Stabilizing Patches
            if resA_name in HYDROPHOBIC_AAS and resB_name in HYDROPHOBIC_AAS:
                hydrophobic_patches.append((resA_info, resB_info))



           # ta = target[a]
           # bb = binder[b]
            #print(
               # f"Target {ta.resname}{ta.resid}:{ta.name}  —  Binder {bb.resname}{bb.resid}:{bb.name}  =  {d[a, b]:.2f} Å")
            # print(f"Target {ta}:  —  Binder {bb}  =  {d[a,b]:.2f} Å")

    sorted_contacts = sorted(list(pairs), key=lambda x: (x[0][0], x[1][0]))
    sorted_patches = sorted(hydrophobic_patches, key=lambda x: (x[0][0], x[1][0]))

    return sorted_contacts, sorted_patches


def count_hbond_like(target_atoms, binder_atoms, cutoff=3.5):
    """Count N/O atom pairs within cutoff. Return (count, unique_pairs_list).
       pairs reported as (target_resid, binder_resid)
    """
    targ_NO = target_atoms.select_atoms("name N O")
    bind_NO = binder_atoms.select_atoms("name N O")
    if len(targ_NO) == 0 or len(bind_NO) == 0:
        return 0, []
    d = distances.distance_array(targ_NO.positions, bind_NO.positions)
    i, j = np.where(d <= cutoff)
    pairs = set()
    for a, b in zip(i, j):
        pairs.add((int(targ_NO[a].resid), int(bind_NO[b].resid)))
    return len(pairs), sorted(pairs)


def count_clashes(universe, cutoff=2.2):
    """Count inter-chain heavy-atom clashes (pairs < cutoff) across all chains.
       Returns number of pairs found.
    """
    # simple approach: consider all atoms, but exclude same-residue short distances by chain separation
    coords = universe.atoms.positions
    if len(coords) == 0:
        return 0
    # We'll do naive O(N^2) by chunking? For small models it is fine.
    from scipy.spatial import cKDTree
    kdt = cKDTree(coords)
    pairs = kdt.query_pairs(r=cutoff)
    # Filter out pairs that are within same residue (same resid & chain) - keep only inter-entity
    clashes = 0
    for a, b in pairs:
        aatom = universe.atoms[a]
        batom = universe.atoms[b]
        # skip hydrogens
        if aatom.element == 'H' or batom.element == 'H' or aatom.name.startswith('H') or batom.name.startswith('H'):
            continue
        # skip same residue
        if (aatom.segid, aatom.resnum, aatom.resname) == (batom.segid, batom.resnum, batom.resname):
            continue
        # NEW: skip immediate neighbors in the same chain
        if aatom.segid == batom.segid:
            if abs(aatom.resid - batom.resid) <= 1:
                continue
        clashes += 1
    return clashes


def extract_sequence(atomgroup):
    """Return one-letter sequence for an MDAnalysis atomgroup's residues (unique residues in order)."""
    # MDAnalysis resid/resname mapping
    three_to_one = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
                    'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
                    'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

    seq = []
    seen = set()
    for res in atomgroup.residues:
        if res.resid in seen:
            continue
        seen.add(res.resid)
        aa = three_to_one.get(res.resname.strip().upper(), 'X')
        seq.append(aa)
    return ''.join(seq)


# ----------------------------- Main design analysis -----------------------------

def analyze_design(pdb_path, target_chain='A', binder_chain='B', add_target_res_offset=0, freesasa_available=False, tmpdir=None):
    """Analyze a single design PDB file. Returns dict of computed metrics."""
    out = {}
    out['pdb'] = pdb_path
    try:
        u = mda.Universe(pdb_path)
    except Exception as e:
        out['error'] = f"MDAnalysis load error: {e}"
        return out

    # Support chain selection: MDAnalysis may use chainIDs as segids depending on file
    # We'll try both selectors (chainID and segid) for robustness
    def select_non_h(target_chain_letter):
        sel1 = u.select_atoms(f"chainID {target_chain_letter} and not name H*")
        if len(sel1) > 0:
            return sel1
        sel2 = u.select_atoms(f"segid {target_chain_letter} and not name H*")
        if len(sel2) > 0:
            return sel2
        # fallback: try chain by residue range? return empty
        return sel1

    target = select_non_h(target_chain)
    binder = select_non_h(binder_chain)

    out['target_atoms'] = len(target)
    out['binder_atoms'] = len(binder)

    # Contacts
   # 3.0A Contacts and Hydrophobic Patches
    pairs3, hypho_list = residue_contacts(target, binder, cutoff=3.0, add_target_res_offset=add_target_res_offset)
    
    # 4.0A Contacts
    pairs4, _ = residue_contacts(target, binder, cutoff=4.0, add_target_res_offset=add_target_res_offset)
    out['n_contacts_3A'] = len(pairs3)
    out['n_contacts_4A'] = len(pairs4)
    out['pairs_3A'] = pairs3
    out['pairs_4A'] = pairs4
    out['hypho'] = hypho_list

    # Unique contacting residues
    target_res_set = set([p[0][0] for p in pairs4])
    binder_res_set = set([p[1][0] for p in pairs4])
    out['n_target_interface_residues'] = len(target_res_set)
    out['n_binder_interface_residues'] = len(binder_res_set)

    # H-bond-like
    hb_count, hb_pairs = count_hbond_like(target, binder, cutoff=3.5)
    out['hbond_like_count'] = hb_count
    out['hbond_pairs'] = hb_pairs

    # Clashes (fast check). This can be slow for huge systems.
    try:
        clashes = count_clashes(u, cutoff=2.2)
    except Exception:
        clashes = None
    out['clash_count'] = clashes


    # Sequences
    out['target_seq'] = extract_sequence(target)
    out['binder_seq'] = extract_sequence(binder)

    return out


# ----------------------------- CLI & main -----------------------------------

def main():
    parser = argparse.ArgumentParser(description='BindCraft analysis pipeline')
    parser.add_argument('--pdb_dir', required=True, help='Directory with design PDBs (Accepted/Ranked)')
    parser.add_argument('--csv', required=True, help='final_design_stats.csv path')
    parser.add_argument('--out_dir', required=True, help='Output directory for analysis')
    parser.add_argument('--target_chain', default='A', help='Target chain letter (default A)')
    parser.add_argument('--binder_chain', default='B', help='Binder chain letter (default B)')
    parser.add_argument('--add_target_res_offset', type=int, default=0, help='Add offset to target residue numbers if PDB numbering starts at higher value')
    parser.add_argument('--use_freesasa', action='store_true', help='Run freesasa to compute dSASA (requires freesasa on PATH)')
    parser.add_argument('--tmpdir', default=None, help='Temporary directory to write freesasa files')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    freesasa_available = False
    if args.use_freesasa and shutil_which('freesasa'):
        freesasa_available = True
    elif args.use_freesasa:
        print('WARNING: freesasa requested but not found on PATH. dSASA will be skipped.')

    # read CSV metrics
    df = pd.read_csv(args.csv)

    # collect PDB files in directory (match design names in CSV if possible)
    pdb_files = sorted([os.path.join(args.pdb_dir, f) for f in os.listdir(args.pdb_dir) if f.lower().endswith(('.pdb', '.cif'))])

    results = []
    for pdb_path in pdb_files:
        print(f'Analyzing {os.path.basename(pdb_path)}...')
        r = analyze_design(pdb_path, target_chain=args.target_chain, binder_chain=args.binder_chain,
                           add_target_res_offset=args.add_target_res_offset,
                           freesasa_available=freesasa_available, tmpdir=args.tmpdir)
        # merge with CSV row (match by design id in filename)
        base = os.path.splitext(os.path.basename(pdb_path))[0]
        base = "_".join(base.split("_")[:-1])
        # try to find matching row in df by 'Design' or by filename substring
        matched_row = None
        if 'Design' in df.columns:
            # some CSVs have Design column like 'b_F_l75_s...' or similar
            m = df[df['Design'].astype(str).str.contains(base)]
            if len(m) == 1:
                matched_row = m.iloc[0].to_dict()
        if matched_row is None:
            # fallback: try to find row whose Design equals base
            m = df[df['Design'].astype(str) == base]
            if len(m) == 1:
                matched_row = m.iloc[0].to_dict()
        if matched_row is None:
            # no match — keep original df row if filename is same as index
            matched_row = {}
        # build final record
        record = {}
        record.update({'design_id': base})
        # include selected CSV columns if present
        for col in ['Average_pLDDT','Average_i_pLDDT','Average_pTM','Average_i_pTM','Average_pAE','Average_i_pAE','Average_dG','Average_dSASA','Average_Binder_pLDDT','Average_n_InterfaceResidues']:
            record[col] = matched_row.get(col, np.nan)
        # include structural outputs
        record['n_contacts_3A'] = r.get('n_contacts_3A')
        record['n_contacts_4A'] = r.get('n_contacts_4A')
        record['n_target_interface_residues'] = r.get('n_target_interface_residues')
        record['n_binder_interface_residues'] = r.get('n_binder_interface_residues')
        record['hbond_like_count'] = r.get('hbond_like_count')
        record['clash_count'] = r.get('clash_count')
        record['dsasa'] = r.get('dsasa')
        record['target_seq'] = r.get('target_seq')
        record['binder_seq'] = r.get('binder_seq')

        results.append(record)

        # write per-design short report
        rpt_lines = []
        rpt_lines.append(f"Design: {base}")
        rpt_lines.append("--- CSV metrics (if available) ---")
        for col in ['Average_pLDDT','Average_i_pLDDT','Average_pTM','Average_i_pAE','Average_i_pTM','Average_pAE','Average_dG','Average_dSASA','Average_Binder_pLDDT','Average_n_InterfaceResidues']:
            rpt_lines.append(f"{col}: {record.get(col)}")
        rpt_lines.append('')
        rpt_lines.append('--- Structural metrics (computed) ---')
        rpt_lines.append(f"n_contacts_3A: {record['n_contacts_3A']}")
        rpt_lines.append(f"n_contacts_4A: {record['n_contacts_4A']}")
        rpt_lines.append(f"n_target_interface_residues: {record['n_target_interface_residues']}")
        rpt_lines.append(f"n_binder_interface_residues: {record['n_binder_interface_residues']}")
        rpt_lines.append(f"hbond_like_count: {record['hbond_like_count']}")
        rpt_lines.append(f"clash_count: {record['clash_count']}")
        rpt_lines.append(f"dsasa: {record['dsasa']}")
        rpt_lines.append('')
        rpt_lines.append('Sequences:')
        rpt_lines.append(f"target_seq: {record['target_seq']}")
        rpt_lines.append(f"binder_seq: {record['binder_seq']}")

        rpt_path = os.path.join(args.out_dir, f"{base}_report.txt")
        with open(rpt_path, 'w') as fh:
            fh.write('\n'.join(rpt_lines))

    # write master CSV
    df_out = pd.DataFrame(results)
    out_csv = os.path.join(args.out_dir, 'bindcraft_analysis_summary.csv')
    df_out.to_csv(out_csv, index=False)
    print(f"Wrote summary CSV: {out_csv}")

    # produce a simple ranking by i_pTM (csv metric) then dsasa then Average_pLDDT
    df_rank = df_out.copy()
    for col in ['Average_i_pTM','dsasa','Average_pLDDT']:
        if col not in df_rank.columns:
            df_rank[col] = np.nan
    df_rank = df_rank.sort_values(by=['Average_i_pTM','dsasa','Average_pLDDT'], ascending=[False, False, False])
    rank_csv = os.path.join(args.out_dir, 'bindcraft_ranked.csv')
    df_rank.to_csv(rank_csv, index=False)
    print(f"Wrote ranked CSV: {rank_csv}")


if __name__ == '__main__':
    main()

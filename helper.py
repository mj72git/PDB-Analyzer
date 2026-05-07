import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import tempfile
import os
import numpy as np
import plotly.express as px
import MDAnalysis as mda
import tempfile, os


############# HELPER FUNCTIONS #############

def extract_contact_residues(pairs):
    """Extract unique residue numbers from a list of contact pairs [( (A,res1),(B,res2) ), ... ]"""
    if not pairs or pairs == "nan":
        return []
    residues = set()
    try:
        for (a_chain, a_res), (b_chain, b_res) in pairs:
            residues.add((a_chain, int(a_res)))
            residues.add((b_chain, int(b_res)))
    except Exception:
        return []
    return list(residues)

def parse_pairs(raw_pairs, target_chain='A', binder_chain='B'):
    """
    Parse BindCraft style contact pairs.
    Example input:
    [((75,'ARG'), (28,'PHE')), ((78,'THR'), (35,'TRP'))]
    Returns:
       [('A',75), ('B',28), ('A',78), ('B',35)]
    """
    import ast

    if raw_pairs is None:
        return []

    # Convert string → Python list
    if isinstance(raw_pairs, str):
        try:
            pairs = ast.literal_eval(raw_pairs)
        except Exception:
            return []
    else:
        pairs = raw_pairs

    residues = set()

    #[((75,'ARG'), (28,'PHE')), ((78,'THR'), (35,'TRP'))]
    for left, right in pairs:
        try:
            # left = (75, 'ARG') → residue number = left[0]
            res_target = int(left[0])
            residues.add((target_chain, res_target))

            # right = (28, 'PHE') → residue number = right[0]
            res_binder = int(right[0])
            residues.add((binder_chain, res_binder))

        except Exception:
            continue

    return list(residues)

import re

def extract_design_id(pdb_filename):

    """
    Works for:
    - b_SLOG_l50_s95426_mpnn9_model2.pdb
    - 5_b_SLOG_l50_s95426_mpnn9_model2.pdb
    - 12_b_SLOG_l50_s95426_mpnn9_model1.pdb
    """
    name = os.path.splitext(os.path.basename(pdb_filename))[0]

    # remove _modelX suffix
    #name = re.sub(r'_model\d+$', '', name)
    name = "_".join(name.split("_")[:-1])

    # remove leading rank (digits + underscore)
    name = re.sub(r'^\d+_', '', name)

    return name


#######################

def highlight_residues(view, residue_list, sphere=False, cartoon_color='red', sphere_radius=1.2):
    """Highlight residues on the 3D structure as cartoon only."""
    for chain, resi in residue_list:
        try:
            view.addStyle({'chain': chain, 'resi': str(resi)}, {'cartoon': {'color': cartoon_color}})
        except Exception:
            continue
#########################

from MDAnalysis.lib.distances import distance_array

def residue_pair_distance(pdb_text, chain1, res1, chain2, res2):
    import MDAnalysis as mda
    import tempfile, os

    # Write pdb to temp file
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as f:
        f.write(pdb_text.encode())
        fname = f.name

    u = mda.Universe(fname)

    try:
        a1 = u.select_atoms(f"chainID {chain1} and resid {res1}")
        a2 = u.select_atoms(f"chainID {chain2} and resid {res2}")
        if len(a1)==0 or len(a2)==0:
            return None
        d = distance_array(a1.positions, a2.positions).min()
        return round(float(d), 2)
    finally:
        os.remove(fname)


########
def format_pairs_with_distance(pairs, pdb_text, target_chain, binder_chain, add_target_res_offset):
    if pairs is None:
        return "None"

    if isinstance(pairs, str):
        try:
            pairs = eval(pairs)
        except:
            return pairs

    out = []
    for (r1, r2) in pairs:
        # r1 = (resid, resname)
        # r2 = (resid, resname)
        res1, name1 = r1
        res1 = res1 - add_target_res_offset
        res2, name2 = r2

        d = residue_pair_distance(
            pdb_text,
            target_chain, res1,
            binder_chain, res2
        )

        if d is None:
            #out.append(f"{target_chain}{res1 + add_target_res_offset}({name1})   –   {binder_chain}{res2}({name2})   :   NA")
            out.append(f"{res1 + add_target_res_offset} ({name1})     –    {res2}({name2})    :    NA")
        else:
            out.append(f" {res1 + add_target_res_offset} ({name1})      ↔      {res2}({name2})    :   {d} Å")

    return "\n".join(out)


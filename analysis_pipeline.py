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


##################################################################################################
################################# Residue contacts #################################################
##################################################################################################
def residue_contacts(target_atoms, binder_atoms, cutoff=4.0, add_target_res_offset=0, shift = 0):
    """
    Returns:
    - contacts: sorted list of ((resA_id, resA_name), (resB_id, resB_name))
    - hydrophobic_patches: list of contacts specifically between two hydrophobic residues
    """
    # 1. Define Hydrophobic Lookup Table

    if len(target_atoms) == 0 or len(binder_atoms) == 0:
        return [], []

    # 2. Calculate Distances
    d = distances.distance_array(target_atoms.positions, binder_atoms.positions)
    i, j = np.where(d <= cutoff)
    
    pairs = set()
    hydrophobic_patches = []
    salt_bridge = []

    # 3. Analyze Contacts
    flag = ' (A) '
    for a, b in zip(i, j):
        # Extract residue info
        if target_atoms[a].resid < shift :
            flag = ' (A) '
            resA_id = int(target_atoms[a].resid + add_target_res_offset)
            resA_name = target_atoms[a].resname
        else:
            flag = ' (B) '
            resA_id = int(target_atoms[a].resid + add_target_res_offset - shift)
            resA_name = target_atoms[a].resname
        
        resB_id = int(binder_atoms[b].resid)
        resB_name = binder_atoms[b].resname
        
        if shift != 0:
            resA_info = (flag, resA_id, resA_name)
            resB_info = (resB_id, resB_name)
        else:
            flag = ' (A) '
            resA_info = (flag, resA_id, resA_name)
            resB_info = (resB_id, resB_name)

        
        # Add to unique set of pairs
        if (resA_info, resB_info) not in pairs:
            pairs.add((resA_info, resB_info))

    sorted_contacts = sorted(list(pairs), key=lambda x: (x[0][0], x[0][1], x[1][0]))
    return sorted_contacts

##################################################################################################
################################# Analyze interactions #################################################
def analyze_interactions(target_atoms, binder_atoms, add_target_res_offset=0,
                         salt_cutoff=4.0, hydrophobic_cutoff=4.0, shift = 0):
    if len(target_atoms) == 0 or len(binder_atoms) == 0:
        return [], [], []

    # Amino Acid Definitions
    POSITIVE_AAS = {'ARG', 'LYS'}
    NEGATIVE_AAS = {'ASP', 'GLU'}

    # Specific sidechain atom names for salt bridges
    CHARGED_SIDECHAIN_ATOMS = {
        'ASP': ['OD1', 'OD2'],
        'GLU': ['OE1', 'OE2'],
        'LYS': ['NZ'],
        'ARG': ['NH1', 'NH2', 'NE']
    }
    hydrophobic_patches = set()
    salt_bridges = set()

##### -------------------------------------------------------------#################################################
##### 1. SALT BRIDGES (Checked specifically between charged sidechain atoms)#################################################
##### -------------------------------------------------------------#################################################
###### Filter selection to charged sidechains only
    target_charged = target_atoms.select_atoms(
        "resname ASP GLU LYS ARG HIS and name OD1 OD2 OE1 OE2 NZ NH1 NH2 NE ND1 NE2")
    binder_charged = binder_atoms.select_atoms(
        "resname ASP GLU LYS ARG HIS and name OD1 OD2 OE1 OE2 NZ NH1 NH2 NE ND1 NE2")

    if target_charged and binder_charged:
        d_salt = distances.distance_array(target_charged.positions, binder_charged.positions)
        i_s, j_s = np.where(d_salt <= salt_cutoff)

        flag = ' (A) '
        for a, b in zip(i_s, j_s):
            atomA = target_charged[a]
            atomB = binder_charged[b]

            resA_name, resB_name = atomA.resname, atomB.resname

            # Verify opposite charges (Positive-Negative)
            is_A_pos_B_neg = (resA_name in POSITIVE_AAS) and (resB_name in NEGATIVE_AAS)
            is_A_neg_B_pos = (resA_name in NEGATIVE_AAS) and (resB_name in POSITIVE_AAS)

            if is_A_pos_B_neg or is_A_neg_B_pos:
                if atomA.resid < shift :
                    flag = ' (A) '
                    resA_info = (flag, (atomA.resid + add_target_res_offset), resA_name, atomA.name)
                    resB_info = ((atomB.resid), resB_name, atomB.name)
                    salt_bridges.add((resA_info, resB_info))
                else:
                    if shift ==0 :
                        flag = '  '
                    else:
                        flag = ' (B) '
                    resA_info = (flag, (atomA.resid + add_target_res_offset - shift), resA_name, atomA.name)
                    resB_info = ((atomB.resid), resB_name, atomB.name)
                    salt_bridges.add((resA_info, resB_info))

####### -------------------------------------------------------------#################################################
####### 2. HYDROPHOBIC CONTACTS (Checked specifically between sidechain non-polar carbons)#################################################
####### -------------------------------------------------------------
######Exclude backbone atoms (N, C, O, CA) to only measure true hydrophobic sidechain contacts #################################################
    
    HYDROPHOBIC_AAS = "VAL ILE LEU PHE MET TRP ALA PRO"
    HYDROPHOBIC_ATOMS = "CB CG* CD* CE* CZ* CH* SD"

    target_hydro = target_atoms.select_atoms(
        f"resname {HYDROPHOBIC_AAS} and name {HYDROPHOBIC_ATOMS}"
    )
    binder_hydro = binder_atoms.select_atoms(
        f"resname {HYDROPHOBIC_AAS} and name {HYDROPHOBIC_ATOMS}"
    )
    
    if target_hydro and binder_hydro:
        d_hydro = distances.distance_array(target_hydro.positions, binder_hydro.positions)
        i_h, j_h = np.where(d_hydro <= hydrophobic_cutoff)

        flag = ' (A) '
        for a, b in zip(i_h, j_h):
            atomA = target_hydro[a]
            atomB = binder_hydro[b]

            if atomA.resid < shift:
                flag = ' (A) '
                resA_info = (flag, (atomA.resid + add_target_res_offset), atomA.resname, atomA.name)
                resB_info = ((atomB.resid), atomB.resname, atomB.name)
                hydrophobic_patches.add((resA_info, resB_info))
            else:
                if shift == 0:
                    flag = '  '
                else:
                    flag = ' (B) '
                resA_info = (flag, (atomA.resid + add_target_res_offset - shift), atomA.resname, atomA.name)
                resB_info = ((atomB.resid), atomB.resname, atomB.name)
                hydrophobic_patches.add((resA_info, resB_info))

    # Sort output
    sorted_patches = sorted(list(hydrophobic_patches), key=lambda x: (x[0][0], x[0][1], x[1][0]))
    sorted_salt_bridges = sorted(list(salt_bridges), key=lambda x: (x[0][0], x[0][1], x[1][0]))

    return sorted_patches, sorted_salt_bridges

##################################################################################################
################################# pi_pi & cation_pi intercations #################################################

def get_ring_geometry(res):

    AROMATIC_RINGS = {
        'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TRP': ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
        'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2']
    }

    if res.resname not in AROMATIC_RINGS:
        return None

    ring_atoms = res.atoms.select_atoms(
        "name " + " ".join(AROMATIC_RINGS[res.resname])
    )

    expected_atoms = len(AROMATIC_RINGS[res.resname])

    # Require all expected ring atoms
    if len(ring_atoms) != expected_atoms:
        return None

    # Coordinates of all ring atoms
    pos = ring_atoms.positions.astype(float)

    # Geometric centroid
    center = pos.mean(axis=0)

    # Centre the coordinates around the centroid
    centered_pos = pos - center

    # Fit a plane using SVD
    _, _, vh = np.linalg.svd(centered_pos)

    # The last right-singular vector is perpendicular to the fitted plane
    normal = vh[-1]

    # Normalise the vector
    normal = normal / np.linalg.norm(normal)

    return {
        'resid': res.resid,
        'resname': res.resname,
        'center': center,
        'normal': normal
    }


def get_cation_centroids(atom_group):
    """Extracts true positive charge centroids for LYS (NZ) and ARG (NE+NH1+NH2)."""
    cations = []
    for res in atom_group.residues:
        if res.resname == 'LYS':
            nz = res.atoms.select_atoms("name NZ")
            if len(nz) > 0:
                cations.append({
                    'resid': res.resid,
                    'resname': 'LYS',
                    'center': nz.positions[0]
                })
        elif res.resname == 'ARG':
            guanidinium = res.atoms.select_atoms("name NE NH1 NH2")
            if len(guanidinium) == 3:
                # Centroid of the entire guanidinium group
                center = guanidinium.positions.mean(axis=0)
                cations.append({
                    'resid': res.resid,
                    'resname': 'ARG',
                    'center': center
                })
    return cations


def analyze_pi_interactions_strict(target_atoms, binder_atoms, pi_dist_cutoff=5.5, cation_dist_cutoff=6.0,
                                   max_offset=1.39, shift = 0, add_target_res_offset = 0):
    AROMATIC_RINGS = {
        'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TRP': ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
        'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2']
    }

    pi_pi_contacts = []
    cation_pi_contacts = []

    target_rings = [get_ring_geometry(r) for r in target_atoms.residues if get_ring_geometry(r)]
    binder_rings = [get_ring_geometry(r) for r in binder_atoms.residues if get_ring_geometry(r)]

    # 1. RIGOROUS PI-PI INTERACTIONS
    for rA in target_rings:
        for rB in binder_rings:
            dist = np.linalg.norm(rA['center'] - rB['center'])
            if dist <= pi_dist_cutoff:
                # Plane normal angle
                dot_prod = np.clip(abs(np.dot(rA['normal'], rB['normal'])), 0.0, 1.0)
                angle_deg = np.degrees(np.arccos(dot_prod))

                # Calculate lateral offset (projection onto ring plane)
                vec_AB = rB['center'] - rA['center']
                proj_dist = abs(np.dot(vec_AB, rA['normal']))  # Height above plane
                offset = np.sqrt(max(0, dist ** 2 - proj_dist ** 2))  # Lateral shift

                if offset <= max_offset:

                    geometry = None

                    if angle_deg <= 30.0:
                        geometry = "Parallel"

                    elif 60.0 <= angle_deg <= 90.0:
                        geometry = "T-shaped (Edge-to-Face)"

                   
                    # Drop intermediate/unclassified angles to prevent false positives
                    if geometry:
                        if shift == 0:
                            pi_pi_contacts.append({
                                'target': f"{rA['resid'] + add_target_res_offset} {rA['resname']}",
                                'binder': f"{rB['resid']} {rB['resname']}",
                                'distance_A': round(dist, 2),
                                'offset_A': round(offset, 2),
                                'angle_deg': round(angle_deg, 1),
                                'geometry': geometry
                            })
                        else:
                            if rA['resid'] < shift:
                                pi_pi_contacts.append({
                                    'target': f" (A) {rA['resid'] + add_target_res_offset} {rA['resname']}",
                                    'binder': f"{rB['resid']} {rB['resname']}",
                                    'distance_A': round(dist, 2),
                                    'offset_A': round(offset, 2),
                                    'angle_deg': round(angle_deg, 1),
                                    'geometry': geometry
                                })
                            else:
                                pi_pi_contacts.append({
                                    'target': f" (B) {rA['resid'] - shift + add_target_res_offset} {rA['resname']}",
                                    'binder': f"{rB['resid']} {rB['resname']}",
                                    'distance_A': round(dist, 2),
                                    'offset_A': round(offset, 2),
                                    'angle_deg': round(angle_deg, 1),
                                    'geometry': geometry
                                })

    # 2. RIGOROUS CATION-PI INTERACTIONS
    target_cations = get_cation_centroids(target_atoms)
    binder_cations = get_cation_centroids(binder_atoms)

    def check_cation_pi_target_cation(cations, rings):
        results = []
        for cat in cations:
            for ring in rings:
                dist = np.linalg.norm(cat['center'] - ring['center'])
                if dist <= cation_dist_cutoff:
                    cat_vec = cat['center'] - ring['center']
                    proj_dist = abs(np.dot(cat_vec, ring['normal']))  # Vertical height
                    offset = np.sqrt(max(0, dist ** 2 - proj_dist ** 2))  # Lateral offset

                    cat_vec_norm = cat_vec / np.linalg.norm(cat_vec)
                    dot_prod = np.clip(abs(np.dot(ring['normal'], cat_vec_norm)), 0.0, 1.0)
                    theta_deg = np.degrees(np.arccos(dot_prod))

                    ######### adding shift & add_target_res_offset ########################################
                    ################################################################################

                    # Must sit above ring face (theta <= 45°) AND have lateral offset <= 2.0 Å
                    if theta_deg <= 45.0 and offset <= max_offset:
                        if shift == 0:
                            results.append({
                                'cation': f"{cat['resid'] + add_target_res_offset} {cat['resname']}",
                                'aromatic': f"{ring['resid']} {ring['resname']}",
                                'distance_A': round(dist, 2),
                                'offset_A': round(offset, 2),
                                'theta_deg': round(theta_deg, 1)
                            })

                        elif cat['resid'] < shift:
                            results.append({
                                'cation': f" (A) {cat['resid'] + add_target_res_offset} {cat['resname']}",
                                'aromatic': f"{ring['resid']} {ring['resname']}",
                                'distance_A': round(dist, 2),
                                'offset_A': round(offset, 2),
                                'theta_deg': round(theta_deg, 1)
                            })

                        else:
                            results.append({
                                'cation': f" (B) {cat['resid'] - shift + add_target_res_offset} {cat['resname']}",
                                'aromatic': f"{ring['resid']} {ring['resname']}",
                                'distance_A': round(dist, 2),
                                'offset_A': round(offset, 2),
                                'theta_deg': round(theta_deg, 1)
                            })

        return results

    ##########################
    def check_cation_pi_binder_cation(cations, rings):
        results = []
        for cat in cations:
            for ring in rings:
                dist = np.linalg.norm(cat['center'] - ring['center'])
                if dist <= cation_dist_cutoff:
                    cat_vec = cat['center'] - ring['center']
                    proj_dist = abs(np.dot(cat_vec, ring['normal']))  # Vertical height
                    offset = np.sqrt(max(0, dist ** 2 - proj_dist ** 2))  # Lateral offset

                    cat_vec_norm = cat_vec / np.linalg.norm(cat_vec)
                    dot_prod = np.clip(abs(np.dot(ring['normal'], cat_vec_norm)), 0.0, 1.0)
                    theta_deg = np.degrees(np.arccos(dot_prod))

                    ######### adding shift & add_target_res_offset ########################################
                    ################################################################################

                    # Must sit above ring face (theta <= 45°) AND have lateral offset <= 2.0 Å
                    if theta_deg <= 45.0 and offset <= max_offset:
                        if shift == 0:
                            results.append({
                                'cation': f"{cat['resid']} {cat['resname']}",
                                'aromatic': f"{ring['resid'] + add_target_res_offset} {ring['resname']}",
                                'distance_A': round(dist, 2),
                                'offset_A': round(offset, 2),
                                'theta_deg': round(theta_deg, 1)
                            })

                        elif ring['resid'] < shift:
                            results.append({
                                'cation': f"{cat['resid']} {cat['resname']}",
                                'aromatic': f" (A) {ring['resid'] + add_target_res_offset} {ring['resname']}",
                                'distance_A': round(dist, 2),
                                'offset_A': round(offset, 2),
                                'theta_deg': round(theta_deg, 1)
                            })
                        else:
                            results.append({
                                'cation': f"{cat['resid']} {cat['resname']}",
                                'aromatic': f" (B) {ring['resid'] - shift + add_target_res_offset} {ring['resname']}",
                                'distance_A': round(dist, 2),
                                'offset_A': round(offset, 2),
                                'theta_deg': round(theta_deg, 1)
                            })

        return results


    cation_pi_contacts.extend(check_cation_pi_target_cation(target_cations, binder_rings))
    cation_pi_contacts.extend(check_cation_pi_binder_cation(binder_cations, target_rings))

    return pi_pi_contacts, cation_pi_contacts

##################################################################################################
################################# interface_contacts_res_numbers  #################################################

def interface_contacts_res_numbers(pairs):
    if len(pairs) == 0:
        return []
    else:
        interface_contacts = []
        for d in pairs:
            interface_contacts.append(d[0][1])
        interface_contacts = sorted(list(set(interface_contacts)))
        return interface_contacts

##################################################################################################
################################# Extract seq three to one  #################################################

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
    return ''.join(seq), len(seq)

##################################################################################################
################################# Binder Analysis  #################################################

def binder_analysis(binder):
    hydrophobic_res_tot = ['ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP']
    positive_res_tot = ['ARG','LYS']
    negative_res_tot = ['ASP','GLU']
    aromatic_res_tot = ['PHE','TYR','TRP','HIS']
    hydrophobic_res = [(res.resid,res.resname) for res in binder.residues if res.resname in hydrophobic_res_tot ]
    positive_res = [(res.resid,res.resname) for res in binder.residues if res.resname in positive_res_tot ]
    negative_res = [(res.resid,res.resname) for res in binder.residues if res.resname in negative_res_tot]
    aromatic_res = [(res.resid,res.resname) for res in binder.residues if res.resname in aromatic_res_tot]


    return hydrophobic_res, positive_res, negative_res, aromatic_res

################################################################################################################
################################################################################################################
# ----------------------------- Main design analysis -----------------------------
################################################################################################################
################################################################################################################

def analyze_design(pdb_path, target_chain='A', binder_chain='B', add_target_res_offset=0, freesasa_available=False, shift = 0,
                   tmpdir=None, dist_threshod=4.0):
    """Analyze a single design PDB file. Returns dict of computed metrics."""
    out = {}
    out['pdb'] = pdb_path
    try:
        u = mda.Universe(pdb_path)
    except Exception as e:
        out['error'] = f"MDAnalysis load error: {e}"
        return out

    # Support chain selection: MDAnalysis may use chainIDs as segids depending on file
    def select_non_h(target_chain_letter):
        sel1 = u.select_atoms(f"chainID {target_chain_letter} and not name H*")
        if len(sel1) > 0:
            return sel1
        sel2 = u.select_atoms(f"segid {target_chain_letter} and not name H*")
        if len(sel2) > 0:
            return sel2
        return sel1


    def pairs_multi(pairs):
        pairs_multi = []
        for i in range(len(pairs)):
            pairs_multi.append((pairs[i][0][0], pairs[i][0][1]))
        pairs_multi = sorted(list(set(pairs_multi)), key=lambda x: (x[0], x[1]))
        return pairs_multi



    target = select_non_h(target_chain)
    binder = select_non_h(binder_chain)

    out['target_atoms'] = len(target)
    out['binder_atoms'] = len(binder)

    # 1. Contacts Calculation (MUST BE DONE FIRST)
    # 3.0A Contacts and Hydrophobic Patches
    pairs3 = residue_contacts(target, binder, cutoff=3.0, add_target_res_offset=add_target_res_offset, shift = shift)
    if shift !=0:
        pairs3_multi = pairs_multi(pairs3)
    else:
        pairs3_multi = []

    # 4.0A Contacts
    pairs4 = residue_contacts(target, binder, cutoff=4.0, add_target_res_offset=add_target_res_offset, shift = shift)
    if shift !=0:
        pairs4_multi = pairs_multi(pairs4)
    else:
        pairs4_multi = []


    hypho_list, salt_bridge =  analyze_interactions(target, binder, add_target_res_offset=add_target_res_offset,
                             salt_cutoff=dist_threshod, hydrophobic_cutoff=dist_threshod, shift = shift)

    if shift != 0:
        hypho_list_multi = pairs_multi(hypho_list)
        salt_bridge_multi = pairs_multi(salt_bridge)
    else:
        hypho_list_multi = []
        salt_bridge_multi = []

    interface_contacts_3A = interface_contacts_res_numbers(pairs3)
    interface_contacts_4A = interface_contacts_res_numbers(pairs4)
    interface_contacts_hypho = interface_contacts_res_numbers(hypho_list)
    interface_contacts_salt_bridge = interface_contacts_res_numbers(salt_bridge)

    pi_pi_contacts, cation_pi_contacts = analyze_pi_interactions_strict(target, binder, pi_dist_cutoff=5.0, cation_dist_cutoff=6.0,
                                   max_offset=1.39, shift = shift, add_target_res_offset = add_target_res_offset)

################################# Continue Binder Analysis  #################################################
    hydrophobic_res, positive_res, negative_res, aromatic_res = binder_analysis(binder)

    hydrophobic_res_fraction = 100 * len(hydrophobic_res) / len(binder.residues)
    positive_res_fraction = 100 * len(positive_res) / len(binder.residues)
    negative_res_fraction = 100 * len(negative_res) / len(binder.residues)
    aromatic_res_fraction = 100 * len(aromatic_res) / len(binder.residues)
    net_charged_estimated = len(positive_res) - len(negative_res)


    hydrophobic_res_interface = []
    for i in range(len(pairs4)):
        if pairs4[i][1] in hydrophobic_res:
            hydrophobic_res_interface.append(pairs4[i][1])

    positive_res_interface = []
    for i in range(len(pairs4)):
        if pairs4[i][1] in positive_res:
            positive_res_interface.append(pairs4[i][1])

    negative_res_interface = []
    for i in range(len(pairs4)):
        if pairs4[i][1] in negative_res:
            negative_res_interface.append(pairs4[i][1])

    aromatic_res_interface = []
    for i in range(len(pairs4)):
        if pairs4[i][1] in aromatic_res:
            aromatic_res_interface.append(pairs4[i][1])

    hydrophobic_res_fraction_interface = 0 if len(pairs4) == 0  else 100 * len(hydrophobic_res_interface) / len(pairs4)
    positive_res_fraction_interface = 0 if len(pairs4) == 0  else 100 * len(positive_res_interface) / len(pairs4)
    negative_res_fraction_interface = 0 if len(pairs4) == 0  else 100 * len(negative_res_interface) / len(pairs4)
    aromatic_res_fraction_interface = 0 if len(pairs4) == 0  else 100 * len(aromatic_res_interface) / len(pairs4)

##################################################################################################
################################# H Bond  #################################################

    def evaluate_pairs_bindcraft(
            donors, acceptors,
            distance_cutoff=3.5,
            angle_cutoff=90.0,
            add_target_res_offset=add_target_res_offset,
            flag='t',
            shift=0):

        ACCEPTOR_PARENTS = {
            'O': 'C',
            'OD1': 'CG',
            'OD2': 'CG',
            'OE1': 'CD',
            'OE2': 'CD',
            'OG': 'CB',
            'OG1': 'CB',
            'OH': 'CZ',
            'ND1': 'CG',
            'NE2': 'CD2',
        }

        results = []

        t = add_target_res_offset if flag == 't' else 0
        b = add_target_res_offset if flag == 'b' else 0

        shift_d = shift if flag == 't' else 0
        shift_a = shift if flag == 'b' else 0

        for d in donors:
            for a in acceptors:

                # 1. Calculate donor-acceptor distance
                dist = np.linalg.norm(d.position - a.position)

                if dist > distance_cutoff:
                    continue

                # 2. Find acceptor parent atom
                parent_name = ACCEPTOR_PARENTS.get(a.name)

                if parent_name is None:
                    continue

                parent_atoms = a.residue.atoms.select_atoms(
                    f'name {parent_name}'
                )

                # Do not accept interactions without the required geometry
                if len(parent_atoms) != 1:
                    continue

                parent_pos = parent_atoms.positions[0]

                # 3. Calculate D-A-Parent angle
                v1 = d.position - a.position
                v2 = parent_pos - a.position

                norm_product = np.linalg.norm(v1) * np.linalg.norm(v2)

                if norm_product == 0:
                    continue

                cosine = np.dot(v1, v2) / norm_product

                angle_deg = np.degrees(
                    np.arccos(np.clip(cosine, -1.0, 1.0))
                )

                # 4. Apply angular criterion
                if angle_deg < angle_cutoff:
                    continue

                # 5. Preserve your original residue numbering logic
                if flag == 't':
                    is_first = d.resid < shift
                else:
                    is_first = a.resid < shift

                if is_first:
                    chain_label = '  ' if shift == 0 else ' (A) '
                    donor_resid = d.resid + t
                    acceptor_resid = a.resid + b
                else:
                    chain_label = '  ' if shift == 0 else ' (B) '
                    donor_resid = d.resid + t - shift_d
                    acceptor_resid = a.resid + b - shift_a

                results.append({
                    'flag_target_chain': chain_label,
                    'donor': (donor_resid, d.resname, d.name),
                    'acceptor': (acceptor_resid, a.resname, a.name),
                    'distance_A': round(float(dist), 2),
                    'angle_deg': round(float(angle_deg), 1),
                })

        return results

    DONOR_NAMES = {
        'N',
        'NE', 'NE1', 'NE2',
        'NH1', 'NH2',
        'ND1', 'ND2',
        'NZ',
        'OG', 'OG1', 'OH',
        'SG',
    }

    ACCEPTOR_NAMES = {
        'O',
        'OD1', 'OD2',
        'OE1', 'OE2',
        'OG', 'OG1', 'OH',
        'ND1', 'NE2',
    }

    # Convert sets to MDAnalysis selection strings
    donor_sel_str = "name " + " ".join(DONOR_NAMES)
    acceptor_sel_str = "name " + " ".join(ACCEPTOR_NAMES)

    # Target Donors & Acceptors
    target_donors = target.select_atoms(donor_sel_str)
    target_acceptors = target.select_atoms(acceptor_sel_str)

    # Binder Donors & Acceptors
    binder_donors = binder.select_atoms(donor_sel_str)
    binder_acceptors = binder.select_atoms(acceptor_sel_str)

    # Direction 1: Target gives H, Binder receives H
    hbonds_target_to_binder = evaluate_pairs_bindcraft(
        donors=target_donors,
        acceptors=binder_acceptors,
        distance_cutoff=3.5,
        angle_cutoff=90.0,add_target_res_offset = add_target_res_offset, flag = 't', shift = shift
    )

    # Direction 2: Binder gives H, Target receives H
    hbonds_binder_to_target = evaluate_pairs_bindcraft(
        donors=binder_donors,
        acceptors=target_acceptors,
        distance_cutoff=3.5,
        angle_cutoff=90.0, add_target_res_offset = add_target_res_offset, flag = 'b', shift = shift
    )

########################################################################################################
########################################################################################################

    out['hbonds_target_to_binder'] = hbonds_target_to_binder
    out['hbonds_binder_to_target'] = hbonds_binder_to_target
    out['Number of All H bonds'] = len(hbonds_target_to_binder) + len(hbonds_binder_to_target)

    #### binder only ###################################
    out['hydrophobic_res'] = hydrophobic_res
    out['hydrophobic_res_fraction'] = round(hydrophobic_res_fraction,2)
    out['hydrophobic_res_interface'] = list(set(hydrophobic_res_interface))
    out['hydrophobic_res_fraction_interface'] = round(hydrophobic_res_fraction_interface,2)

    out['positive_res'] = positive_res
    out['positive_res_fraction'] = round(positive_res_fraction,2)
    out['positive_res_interface'] = list(set(positive_res_interface))
    out['positive_res_fraction_interface'] = round(positive_res_fraction_interface,2)

    out['negative_res'] = negative_res
    out['negative_res_fraction'] = round(negative_res_fraction,2)
    out['negative_res_interface'] = list(set(negative_res_interface))
    out['negative_res_fraction_interface'] = round(negative_res_fraction_interface,2)

    out['aromatic_res'] = aromatic_res
    out['aromatic_res_fraction'] = round(aromatic_res_fraction,2)
    out['aromatic_res_interface'] = list(set(aromatic_res_interface))
    out['aromatic_res_fraction_interface'] = round(aromatic_res_fraction_interface,2)

    out['net_charged_estimated'] = net_charged_estimated
    #########################################
    out['pi_pi_contacts'] = pi_pi_contacts
    out['cation_pi_contacts'] = cation_pi_contacts
    out['n_pi_pi_contacts'] = len(pi_pi_contacts)
    out['n_cation_pi_contacts'] = len(cation_pi_contacts)

    n_hypho_list = len(hypho_list)
    n_salt_bridge = len(salt_bridge)
    out['n_contacts_3A'] = len(pairs3)
    out['n_contacts_4A'] = len(pairs4)
    out['pairs_3A'] = pairs3
    out['pairs3_multi'] = pairs3_multi
    out['pairs_4A'] = pairs4
    out['pairs4_multi'] = pairs4_multi
    out['3A_Residues'] = interface_contacts_3A
    out['4A_Residues'] = interface_contacts_4A
    out['n_contacts_3A_Residues'] = len(interface_contacts_3A)
    out['n_contacts_4A_Residues'] = len(interface_contacts_4A)

    out['hypho'] = hypho_list
    out['hypho_multi'] = hypho_list_multi
    out['salt_bridge'] = salt_bridge
    out['salt_bridge_multi'] = salt_bridge_multi

    out['hypho_Residues'] = interface_contacts_hypho
    out['salt_bridge_Residues'] = interface_contacts_salt_bridge
    out['number_hydrophobic_contacts_Residues'] = len(interface_contacts_hypho)
    out['number_salt_bridge_Residues'] = len(interface_contacts_salt_bridge)
    out['n_hypho_list'] = n_hypho_list
    out['n_salt_bridge'] = n_salt_bridge


    # 2. Unique contacting residues
    target_res_set = set([p[0][0] for p in pairs3])
    binder_res_set = set([p[1][0] for p in pairs3])
    out['n_target_interface_residues'] = len(target_res_set)
    out['n_binder_interface_residues'] = len(binder_res_set)


    # 6. Sequences extraction
    out['target_seq'],_ = extract_sequence(target)
    out['binder_seq'],_ = extract_sequence(binder)
    _,out['binder_length'] = extract_sequence(binder)

    return out

######    FINISHED ##################################

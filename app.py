import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import tempfile
import os
import numpy as np
import plotly.express as px

###### MODIFIED #####
if "analysis_done" not in st.session_state:
    st.session_state.analysis_done = False
####################
try:
    import py3Dmol
    py3dmol_available = True
except ImportError:
    py3dmol_available = False

from analysis_pipeline import analyze_design
from helper import *
st.set_page_config(page_title="PDB Analyser", layout="wide")
st.image('uoa.jpg', width=200)
st.title("PDB Analyser Web Application")
################################################################################
################################################################################
def format_pairs(pairs):
    if not pairs or pairs == 'nan':
        return "None"
    try:
        #return " ".join([f"{a[0]} {a[1]} ↔ {b[0]} {b[1]}" for a, b in pairs])
        #return " ".join([f"{a[0]} {a[1]} {a[2]}   ↔   {b[0]} {b[1]} {b[2]} \n" for a, b in pairs])
        return " ".join([f"{a[0]} {a[1]}  ({a[2]})    ↔   {b[0]} {b[1]}   ({b[2]})  \n" for a, b in pairs])

    except Exception:
        return str(pairs)

def format_pairs_contacts_only(pairs):
    if not pairs or pairs == 'nan':
        return "None"
    try:
        return " ".join([f"{a[0]} {a[1]} ↔ {b[0]} {b[1]} \n" for a, b in pairs])
        #return " ".join([f"{a[0]} {a[1]} ** {a[2]}    ↔   {b[0]} {b[1]}  ** {b[2]}  \n" for a, b in pairs])

    except Exception:
        return str(pairs)

def format_pairs_h_bonds_b_to_t(pairs):
    if not pairs or pairs == 'nan':
        return "None"
    try:
        return " ".join([f" {a['donor'][0]} {a['donor'][1]} ({a['donor'][2]})    ->    {a['acceptor'][0]} {a['acceptor'][1]} ({a['acceptor'][2]})    |    distance = {a['distance_A']} Å    |    Angle = {a['angle_deg']}  degree  \n" for a in pairs])
    except Exception:
        return str(pairs)

def format_pairs_h_bonds_t_to_b(pairs):
    if not pairs or pairs == 'nan':
        return "None"
    try:
        return " ".join([f" {a['donor'][0]} {a['donor'][1]} ({a['donor'][2]})    ->    {a['acceptor'][0]} {a['acceptor'][1]} ({a['acceptor'][2]})    |    distance = {a['distance_A']} Å     |    Angle = {a['angle_deg']}  degree  \n" for a in pairs])
    except Exception:
        return str(pairs)


def format_pairs_pi_pi(pairs):
    if not pairs or pairs == 'nan':
        return "None"
    try:
        return " ".join([f" (target)  {a['target']}   ->     (binder)  {a['binder']}  |  distance = {round(a['distance_A'],2)} Å  |  offset = {a['offset_A']} | angle = {a['angle_deg']} degree  |  geometry = {a['geometry']}  \n" for a in pairs])
    except Exception:
        return str(pairs)

def format_pairs_cation_pi(pairs):
    if not pairs or pairs == 'nan':
        return "None"
    try:
        return " ".join([f" Cation : {a['cation']}   ->   Aromatic : {a['aromatic']}  |  distance = {round(a['distance_A'],2)} Å  |  offset = {a['offset_A']} | angle = {a['theta_deg']} degree \n" for a in pairs])
    except Exception:
        return str(pairs)


def format_pairs_binder(pairs):
    if not pairs or pairs == 'nan':
        return "None"
    try:
        return " ".join([f" {a} \n" for a in pairs])
    except Exception:
        return str(pairs)


################################################################################
################################################################################

# session_state defaults
if 'df_out' not in st.session_state:
    st.session_state.df_out = None
if 'df_rank' not in st.session_state:
    st.session_state.df_rank = None
if "df_final" not in st.session_state:
  st.session_state.df_final = None
if 'pdb_map' not in st.session_state:
    st.session_state.pdb_map = {}


################################### SIDEBAR SETTINGS ##########################

st.sidebar.header("Settings")
target_chain = st.sidebar.text_input("Target Chain Letter", value="A")
binder_chain = st.sidebar.text_input("Binder Chain Letter", value="B")
#dist_threshod = st.sidebar.number_input("Salt bridge/Hydrophobic contacts distance cut off (Å) ", value=4)
check_box_family = st.sidebar.checkbox('In all designs, target protein is same', value = True)
#max_dist_threshod = st.sidebar.number_input("Maximum atomic distance ", value=4)

add_target_res_offset = st.sidebar.number_input("Target Residue Offset ", value=0)
st.sidebar.markdown("##### (e.g. your Output  target chain starts from Residue 1 but your initial target chain, starts from 24. you should enter 23)")
blankk = st.sidebar.header("")

file_path = "help.txt"
with open(file_path, "r") as ff:
    txt_content = ff.read()
st.sidebar.download_button(
    label="Download HELP",
    data=txt_content,
    file_name="help.txt",
    mime="text/plain"
)
mj = st.sidebar.header("App created by MJ Shadfar")
st.sidebar.markdown("#### currently a PhD candidate in A/Prof Jane Allison’s group, School of Biological Science, University of Auckland. ")
st.sidebar.caption(" [Linkedin Profile](https://www.linkedin.com/in/mj-shadfar-3919b5b8/)")
st.sidebar.write(" [Github Repository](https://github.com/mj72git/PDB-Analyzer)")
st.sidebar.caption("PDB Analyser v1.3.1")
st.sidebar.caption("")
st.sidebar.caption("The last modify : 22 Sep 2026")
# st.table(
#     {
#         ":material/folder: Project": "**Streamlit** - The fastest way to build data apps",
#         ":material/code: Repository": "[github.com/streamlit/streamlit](https://github.com/streamlit/streamlit)",
#         ":material/new_releases: Version": ":gray-badge[1.45.0]",
#         ":material/license: License": ":green-badge[Apache 2.0]",
#         ":material/group: Maintainers": ":blue-badge[Core Team] :violet-badge[Community]",
#     },
#     border="horizontal")

############################# FILE UPLOAD & ANALYSIS ########################
if not st.session_state.analysis_done:
    st.subheader("")
    st.subheader("This App is created for analysing pdb files generated from Protein Design software like AlphaFold or BindCraft. (Only pdb files with target & binder) ")
    #st.subheader("")
    st.image('image.jpg', width = 500)
    st.subheader("")
    st.subheader("Upload Files")
    pdb_files = st.file_uploader("Upload PDB files", type=["pdb"], accept_multiple_files=True)
    #csv_file = st.file_uploader("Upload design metrics CSV. (final_design_stats.csv)", type="csv")

    #if pdb_files and csv_file:
    if pdb_files:
        st.success("Files uploaded successfully.")
        if st.button("Run Analysis"):
            tmpdir = tempfile.mkdtemp()
            #freesasa_available = (os.system("which freesasa > /dev/null") == 0)
            #csv_path = os.path.join(tmpdir, "metrics.csv")
            #with open(csv_path, "wb") as f:
                #f.write(csv_file.read())
            #df_metrics = pd.read_csv(csv_path)

            results = []
            progress = st.progress(0)
            for i, pdb in enumerate(pdb_files):
                pdb_path = os.path.join(tmpdir, pdb.name)
                with open(pdb_path, "wb") as f:
                    f.write(pdb.read())

                r = analyze_design(pdb_path, target_chain=target_chain, binder_chain=binder_chain,
                                   add_target_res_offset=add_target_res_offset,tmpdir=tmpdir,dist_threshod=4)

                #base = os.path.splitext(os.path.basename(pdb_path))[0]  #omitt the .pdb
                #base = "_".join(base.split("_")[:-1]) #or os.path.splitext(os.path.basename(pdb_path))[0]


                ######## changed!
                #base = extract_design_id(pdb.name)
                # This removes the '.pdb' extension but leaves the rest of the filename completely intact
                base = os.path.splitext(pdb.name)[0]
                ########

                
                #base = "_".join(base.split("_")[:-1]) or os.path.splitext(os.path.basename(pdb_path))[0]
                #matched_row = {}

                #if 'Design' in df_metrics.columns:

                    #m = df_metrics[df_metrics['Design'].astype(str).str.contains(base)]
                    #if len(m) == 1:
                        #matched_row = m.iloc[0].to_dict()
                    #else:
                        #m = df_metrics[df_metrics['Design'].astype(str) == base]
                        #if len(m) == 1:
                           # matched_row = m.iloc[0].to_dict()

                try:
                    with open(pdb_path, "r") as fh:
                        pdb_text = fh.read()
                except UnicodeDecodeError:
                    with open(pdb_path, "r", encoding="latin-1") as fh:
                        pdb_text = fh.read()
                st.session_state.pdb_map[base] = pdb_text

                record = {'design_id': base}
                #for col in ['Average_pLDDT','Average_i_pLDDT','Average_pTM','Average_i_pTM','Average_pAE','Average_i_pAE','Average_dG','Average_dSASA','Average_Binder_pLDDT','Average_n_InterfaceResidues']:
                    #record[col] = matched_row.get(col, np.nan)
                record.update({
                    'binder_length' : r.get('binder_length'),
                    'Average_pLDDT' : r.get('Average-pLDDT'),
                    #'Average_i_pLDDT' : r.get('Average-i-pLDDT'),
                    #'Average_binder_pLDDT' : r.get('Average-binder-pLDDT'),
                    'n_contacts_3A_Atom': r.get('n_contacts_3A'),
                    'n_contacts_4A_Atom': r.get('n_contacts_4A'),
                    #'n_target_interface_residues': r.get('n_target_interface_residues'),
                    #'n_binder_interface_residues': r.get('n_binder_interface_residues'),
                   # 'hbond_like_count': r.get('hbond_like_count'),
                   # 'clash_count': r.get('clash_count'),
                   # 'dsasa': r.get('dsasa'),
                    'target_seq': r.get('target_seq'),
                    'pairs_3A': r.get('pairs_3A'),
                    '3A_Residues': r.get('3A_Residues'),
                    'n_contacts_3A_Residues' : r.get('n_contacts_3A_Residues'),
                    'pairs_4A': r.get('pairs_4A'),
                    '4A_Residues': r.get('4A_Residues'),
                    'n_contacts_4A_Residues': r.get('n_contacts_4A_Residues'),
                  #  'hbond_pairs': r.get('hbond_pairs')
                    'hydrophobic_contacts' : r.get('hypho'),
                    'n_hydrophobic_contacts_ATOM': r.get('n_hypho_list'),
                    'hydrophobic_contacts_Residues' : r.get('hypho_Residues'),
                    'n_hydrophobic_contacts_Residues' : r.get('number_hydrophobic_contacts_Residues'),
                    'salt_bridge' : r.get('salt_bridge'),
                    'n_salt_bridge_ATOM': r.get('n_salt_bridge'),
                    'salt_bridge_contacts_Residues': r.get('salt_bridge_Residues'),
                    'n_salt_bridge_Residues' : r.get('number_salt_bridge_Residues'),
                    'pi_pi_contacts' : r.get('pi_pi_contacts'),
                    'n_pi_pi_contacts' : r.get('n_pi_pi_contacts'),
                    'cation_pi_contacts' : r.get('cation_pi_contacts'),
                    'n_cation_pi_contacts' : r.get('n_cation_pi_contacts'),
                    'hbonds_target_to_binder' : r.get('hbonds_target_to_binder'),
                    'hbonds_binder_to_target' : r.get('hbonds_binder_to_target'),
                    'Number of All H bonds' : r.get('Number of All H bonds'),
                    'binder_seq': r.get('binder_seq'),
                    'binder_hydrophobic_res': r.get('hydrophobic_res'),
                    'binder_hydrophobic_res_fraction': r.get('hydrophobic_res_fraction'),
                    'binder_hydrophobic_res_interface': r.get('hydrophobic_res_interface'),
                    'binder_hydrophobic_res_fraction_interface': r.get('hydrophobic_res_fraction_interface'),
                    'binder_positive_res': r.get('positive_res'),
                    'binder_positive_res_fraction': r.get('positive_res_fraction'),
                    'binder_positive_res_interface': r.get('positive_res_interface'),
                    'binder_positive_res_fraction_interface': r.get('positive_res_fraction_interface'),
                    'binder_negative_res': r.get('negative_res'),
                    'binder_negative_res_fraction': r.get('negative_res_fraction'),
                    'binder_negative_res_interface': r.get('negative_res_interface'),
                    'binder_negative_res_fraction_interface': r.get('negative_res_fraction_interface'),
                    'binder_aromatic_res': r.get('aromatic_res'),
                    'binder_aromatic_res_fraction': r.get('aromatic_res_fraction'),
                    'binder_aromatic_res_interface': r.get('aromatic_res_interface'),
                    'binder_aromatic_res_fraction_interface': r.get('aromatic_res_fraction_interface')

                })
                results.append(record)
                progress.progress((i+1)/len(pdb_files))

            st.session_state.df_out = pd.DataFrame(results)
            df_rank = st.session_state.df_out.copy()
            #for col in ['Average_i_pTM','Average_dSASA','Average_pLDDT']:
                #if col not in df_rank.columns:
                    #df_rank[col] = np.nan
            st.session_state.df_rank = df_rank.sort_values(by=['n_contacts_3A_Atom','Average_pLDDT','n_hydrophobic_contacts_ATOM','n_salt_bridge_ATOM','n_contacts_4A_Atom'], ascending=[False, False, False, False, False])
            st.session_state.analysis_done = True
            st.rerun()

############################ AFTER ANALYSIS ####################################

if st.session_state.analysis_done and st.session_state.df_out is not None:
    tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs(["Summary", "Visualizations", "Details", "Binder Analysis",
                                                        "Binder Analysis Overall","Overall", "Filters (Customization)"])
############################################################################################################################
############################################################################################################################
############################################################################################################################
    with tab1:
        st.subheader("Analysis Summary")
        df_display = st.session_state.df_out.copy()
        for col in ['binder_length','binder_hydrophobic_res','binder_hydrophobic_res_fraction','binder_hydrophobic_res_interface','binder_hydrophobic_res_fraction_interface',
                    'binder_positive_res','binder_positive_res_fraction', 'binder_positive_res_interface','binder_positive_res_fraction_interface','binder_negative_res','binder_negative_res_fraction',
                    'binder_negative_res_interface','binder_negative_res_fraction_interface','hbonds_target_to_binder','hbonds_binder_to_target',
                    'binder_aromatic_res','binder_aramotic_res_fraction','binder_aromatic_res_interface', 'binder_aromatic_res_fraction_interface','Average_pLDDT','pairs_3A', 'n_salt_bridge_Residues','n_contacts_3A_Residues','n_contacts_4A_Residue','pairs_4A', 'hydrophobic_contacts', 'salt_bridge', 'pi_pi_contacts', 'cation_pi_contacts']:
            if col in df_display.columns:
                df_display[col] = df_display[col].apply(lambda x: str(x))
        st.dataframe(df_display)

        st.subheader("Ranked Designs")
        df_rank_display = st.session_state.df_rank.copy()
        for col in ['binder_length', 'binder_hydrophobic_res', 'binder_hydrophobic_res_fraction',
                    'binder_hydrophobic_res_interface', 'binder_hydrophobic_res_fraction_interface',
                    'binder_positive_res', 'binder_positive_res_fraction', 'binder_positive_res_interface',
                    'binder_positive_res_fraction_interface', 'binder_negative_res', 'binder_negative_res_fraction',
                    'binder_negative_res_interface', 'binder_negative_res_fraction_interface','hbonds_target_to_binder','hbonds_binder_to_target',
                    'binder_aromatic_res', 'binder_aramotic_res_fraction', 'binder_aromatic_res_interface',
                    'binder_aromatic_res_fraction_interface', 'Average_pLDDT', 'pairs_3A', 'n_salt_bridge_Residues',
                    'n_contacts_3A_Residues', 'n_contacts_4A_Residue', 'pairs_4A', 'hydrophobic_contacts',
                    'salt_bridge', 'pi_pi_contacts', 'cation_pi_contacts']:
            if col in df_rank_display.columns:
                df_rank_display[col] = df_rank_display[col].apply(lambda x: str(x))
        st.dataframe(df_rank_display)

        

      

        st.download_button("Download Summary CSV", st.session_state.df_out.to_csv(index=False), "summary.csv")
        st.download_button("Download Ranked CSV", st.session_state.df_rank.to_csv(index=False), "ranked.csv")

############################################################################################################################
############################################################################################################################
############################################################################################################################
    with tab2:
        st.subheader("Visualizations")

        # ensure df is defined
        df = st.session_state.df_out

        # -------------------------------
        # BAR PLOT (4Å contacts)
        # -------------------------------


        if df is not None and df.shape[0] > 0:
            fig_contacts = px.bar(
                df,
                x='design_id',
                y='n_contacts_3A_Residues',
                title='Number of Contacts (3Å)'
            )
            st.plotly_chart(fig_contacts, use_container_width=False)
        else:
            st.info("No data to plot. Please run analysis first.")

        st.write("------------------------------------------------------------------------------------------")
        if df is not None and df.shape[0] > 0:
            fig_contacts = px.bar(
                df,
                x='design_id',
                y='n_contacts_4A_Residues',
                title='Number of Contacts (4Å)'
            )
            st.plotly_chart(fig_contacts, use_container_width=True)
        else:
            st.info("No data to plot. Please run analysis first.")

        st.write("------------------------------------------------------------------------------------------")
        if df is not None and df.shape[0] > 0:
            fig_contacts = px.bar(
                df,
                x='design_id',
                y='n_hydrophobic_contacts_Residues',
                title='Number of hydrophobic Contacts '
            )
            st.plotly_chart(fig_contacts, use_container_width=True)
        else:
            st.info("No data to plot. Please run analysis first.")

        st.write("------------------------------------------------------------------------------------------")
        if df is not None and df.shape[0] > 0:
            fig_contacts = px.bar(
                df,
                x='design_id',
                y='n_salt_bridge_Residues',
                title='Number of salt bridge contacts '
            )
            st.plotly_chart(fig_contacts, use_container_width=True)
        else:
            st.info("No data to plot. Please run analysis first.")

        st.write("------------------------------------------------------------------------------------------")
        if df is not None and df.shape[0] > 0:
            fig_contacts_pi_pi = px.bar(
                df,
                x='design_id',
                y='n_pi_pi_contacts',
                title='Number of Pi_Pi contacts '
            )
            st.plotly_chart(fig_contacts_pi_pi, use_container_width=True)
        else:
            st.info("No data to plot. Please run analysis first.")

        st.write("------------------------------------------------------------------------------------------")
        if df is not None and df.shape[0] > 0:
            fig_contacts_cation_pi = px.bar(
                df,
                x='design_id',
                y='n_cation_pi_contacts',
                title='Number of Cation_Pi contacts '
            )
            st.plotly_chart(fig_contacts_cation_pi, use_container_width=True)
        else:
            st.info("No data to plot. Please run analysis first.")

        st.write("------------------------------------------------------------------------------------------")
        if df is not None and df.shape[0] > 0:
            fig_contacts_hb = px.bar(
                df,
                x='design_id',
                y='Number of All H bonds',
                title='Number of all H bonds '
            )
            st.plotly_chart(fig_contacts_hb, use_container_width=True)
        else:
            st.info("No data to plot. Please run analysis first.")


        # -------------------------------
        # 3D STRUCTURE VIEWER
        # -------------------------------
        st.subheader("3D Structure Viewer")

        if not py3dmol_available:
            st.error("py3Dmol is not installed.")
        else:
            if df is None or df.shape[0] == 0:
                st.info("No designs available. Please run analysis first.")
            else:
                # Design selection
                selected_design = st.selectbox("Select a design", df['design_id'].tolist())
                row = df[df['design_id'] == selected_design].iloc[0]

                # Handle raw_pairs
                raw_pairs = row.get('pairs_4A')
                if isinstance(raw_pairs, str) and raw_pairs.strip().lower() in ["nan", "", "none"]:
                    raw_pairs = None

                # Parse contacts
                contacts = parse_pairs(raw_pairs, target_chain=target_chain, binder_chain=binder_chain)

                # Apply target offset and normalize chains
                contacts_offset = []
                t_chain = str(target_chain).strip().upper()
                b_chain = str(binder_chain).strip().upper()

                # pairs_4A are ALREADY offset correctly


                for c, r in contacts:
                    try:
                        chain_norm = str(c).strip().upper()
                        resi = int(r)
                        if chain_norm == t_chain:
                            resi -= int(add_target_res_offset)
                        contacts_offset.append((chain_norm, resi))
                    except:
                        continue

                #contacts_offset = [
                    #(str(c).strip().upper(), int(r))
                    #for c, r in contacts
                #]

                # Prepare highlight lists for target and binder
                highlight_target = [(c, r) for c, r in contacts_offset if c == t_chain]
                highlight_binder = [(c, r) for c, r in contacts_offset if c == b_chain]

                pdb_text = st.session_state.pdb_map[selected_design]

                # Chain visibility + colors
                show_target = st.checkbox("Show Target Chain", value=True, key=f"show_target_{selected_design}")
                show_binder = st.checkbox("Show Binder Chain", value=True, key=f"show_binder_{selected_design}")
                highlight_binding = st.checkbox("Highlight binding spot", value=True,
                                                key=f"highlight_{selected_design}")
                if highlight_binding:
                    show_sphere = st.checkbox("Show CPK sphere", value=False, key=f"show_sphere_{selected_design}")
                target_color = st.selectbox("Target Color", ["limegreen", "red", "orange", "magenta", "yellow", "cyan"],
                                            index=0, key=f"tcolor_{selected_design}")
                binder_color = st.selectbox("Binder Color",
                                            ["orange", "red", "deepskyblue", "magenta", "yellow", "cyan"], index=0,
                                            key=f"bcolor_{selected_design}")

                # Create viewer
                view = py3Dmol.view(width=900, height=500)
                view.addModel(pdb_text, 'pdb')
                view.setBackgroundColor('0x303030')

                # Cartoon visibility control
                if show_target:
                    view.setStyle({'chain': t_chain}, {'cartoon': {'color': target_color}})
                else:
                    view.setStyle({'chain': t_chain}, {'cartoon': {'opacity': 0.0}})

                if show_binder:
                    view.setStyle({'chain': b_chain}, {'cartoon': {'color': binder_color}})
                else:
                    view.setStyle({'chain': b_chain}, {'cartoon': {'opacity': 0.0}})

                # Contact highlighting, only for visible chains and if checkbox is checked
                # Contact highlighting: show interface residues as STICKS
                # Contact highlighting: show interface residues as STICKS + keep cartoon
                if highlight_binding:

                    # Target interface
                    if show_target and len(highlight_target) > 0:
                        for chain, resi in highlight_target:
                            view.addStyle(
                                {'chain': chain, 'resi': resi},
                                {'stick': {'colorscheme': 'magentaCarbon', 'radius': 0.25}}
                            )

                    # Binder interface
                    if show_binder and len(highlight_binder) > 0:
                        for chain, resi in highlight_binder:
                            view.addStyle(
                                {'chain': chain, 'resi': resi},
                                {'stick': {'colorscheme': 'yellowCarbon', 'radius': 0.25}}
                            )

                    if show_sphere:
                        target_sphere_color = st.selectbox("Target Sphere Color", ["deepskyblueCarbon", "redCarbon", "orangeCarbon", "magentaCarbon", "yellowCarbon", "cyanCarbon", "spectrum"],
                                            index=2, key=f"tscolor_{selected_design}")
                        binder_sphere_color = st.selectbox("binder Sphere Color", ["limegreenCarbon", "redCarbon", "orangeCarbon", "magentaCarbon", "yellowCarbon", "cyanCarbon", "spectrum"],
                                            index=0, key=f"bscolor_{selected_design}")
                        for chain, resi in highlight_target:
                            view.addStyle(
                                {'chain': chain, 'resi': resi},
                                {'sphere': {'colorscheme': target_sphere_color,'scale':0.8}})
                                #{'sphere': {'colorscheme': 'greenCarbon','scale':0.8}})
                                #{'sphere': {'colorscheme': 'spectrum','scale':0.8}})

                        for chain, resi in highlight_binder:
                            view.addStyle(
                                {'chain': chain, 'resi': resi},
                                {'sphere': {'colorscheme': binder_sphere_color,'scale':0.8}})
                                #{'sphere': {'colorscheme': 'deepskyblueCarbon','scale':0.8}})
                                #{'sphere': {'colorscheme': 'greenCarbon','scale':0.8}})

                view.zoomTo()
                
                components.html(view._make_html(), height=500, width=900)

############################################################################################################################
############################################################################################################################
############################################################################################################################
    with tab3:
        st.subheader("Per-Design Details")
        for _, row in st.session_state.df_out.iterrows():
            with st.expander(f"Details for {row['design_id']}"):
                st.write("**Sequences**")
                st.text(f"Target: {row['target_seq']}")
                st.text(f"Binder: {row['binder_seq']}")
                pdb_text = st.session_state.pdb_map[row['design_id']]
                st.text('')
                #d = format_pairs_with_distance(row['pairs_3A'], pdb_text)
                #st.write("**Contacts (2Å)**")
               # st.text(format_pairs(row['pairs_2A']))
                st.write("**H-bonds**")
                #st.text(f"hbonds_target_to_binder : {row['hbonds_target_to_binder']}")
                st.write("***Target -> Binder***")
                st.text(format_pairs_h_bonds_t_to_b(row['hbonds_target_to_binder']))
                st.text('')
                st.write("***Binder -> Target***")
                st.text(format_pairs_h_bonds_b_to_t(row['hbonds_binder_to_target']))

                #st.text(f"hbonds_binder_to_target : {row['hbonds_binder_to_target']}")
                st.write("------------------------------------------------------------------")
                st.write("**Contacts (3Å)**")
                st.write("***(Target Residues  ↔  binder residues)***")
                st.text(format_pairs_contacts_only(row['pairs_3A']))
                #st.text(format_pairs_with_distance(row['pairs_3A'], pdb_text, target_chain, binder_chain, add_target_res_offset))
                st.write("------------------------------------------------------------------")
                st.write("**hydrophobic contacts**")
                st.text(format_pairs(row['hydrophobic_contacts']))
                st.write("------------------------------------------------------------------")
                st.write("**Salt Bridge**")
                st.text(format_pairs(row['salt_bridge']))
                st.write("------------------------------------------------------------------")
                st.write("**pi_pi_contacts**")
                st.text(format_pairs_pi_pi(row['pi_pi_contacts']))
                #st.text(row['pi_pi_contacts'])
                st.write("------------------------------------------------------------------")
                st.write("**cation_pi_contacts**")
                st.text(format_pairs_cation_pi(row['cation_pi_contacts']))
                #st.text(row['cation_pi_contacts'])
                st.write("------------------------------------------------------------------")
                st.write("------------------------------------------------------------------")
                st.write("**Contacts (4Å)**")
                st.write("***(Target Residues   ↔   binder residues)***")
                st.text(format_pairs_contacts_only(row['pairs_4A']))
                #st.text(format_pairs_with_distance(row['pairs_4A'], pdb_text, target_chain, binder_chain, add_target_res_offset))
                st.write("------------------------------------------------------------------")
############################################################################################################################
############################################################################################################################
############################################################################################################################
    with tab4:
        st.subheader("Per-Binder Details")
        for _, row in st.session_state.df_out.iterrows():
            with st.expander(f"Details for {row['design_id']} -----Binder only"):
                st.write('**Sequence**')
                st.text(f" {row['binder_seq']}")
                st.text(f"length: {row['binder_length']}")
                st.write("------------------------------------------------------------------")
                st.write('**Hydrophobicity/Hydrophilicity Analysis (All aa of the binder):**')
                st.text(f"Hydrophobic: {row['binder_hydrophobic_res_fraction']} %")
                st.text(f"Basic: {row['binder_positive_res_fraction']} %")
                st.text(f"Acidic: {row['binder_negative_res_fraction']} %")
                st.text(f"Aromatic: {row['binder_aromatic_res_fraction']} %")
                st.write("------------------------------------------------------------------")
                st.write('**Hydrophobicity/Hydrophilicity Analysis (Interface Only):**')
                st.text(f"Hydrophobic: {row['binder_hydrophobic_res_fraction_interface']} %")
                st.text(f"Basic: {row['binder_positive_res_fraction_interface']} %")
                st.text(f"Acidic: {row['binder_negative_res_fraction_interface']} %")
                st.text(f"Aromatic: {row['binder_aromatic_res_fraction_interface']} %")
                st.write("------------------------------------------------------------------")
                st.write('**Hydrophobicity Analysis:**')
                st.text(format_pairs_binder(row['binder_hydrophobic_res']))
                st.write("-------------------------------------------------------------------")
                st.write('**Acidic Analysis:**')
                st.text(format_pairs_binder(row['binder_negative_res']))
                st.write("------------------------------------------------------------------")
                st.write('**Basic Analysis:**')
                st.text(format_pairs_binder(row['binder_positive_res']))
                st.write("------------------------------------------------------------------")
                st.write('**Aromatic Analysis:**')
                st.text(format_pairs_binder(row['binder_aromatic_res']))
                st.write("------------------------------------------------------------------")
                st.write("------------------------------------------------------------------")
                st.write('**Hydrophobicity Analysis (Interface):**')
                st.text(format_pairs_binder(row['binder_hydrophobic_res_interface']))
                st.write("------------------------------------------------------------------")
                st.write('**Acidic Analysis (Interface):**')
                st.text(format_pairs_binder(row['binder_negative_res_interface']))
                st.write("------------------------------------------------------------------")
                st.write('**Basic Analysis (Interface):**')
                st.text(format_pairs_binder(row['binder_positive_res_interface']))
                st.write("------------------------------------------------------------------")
                st.write('**Aromatic Analysis (Interface):**')
                st.text(format_pairs_binder(row['binder_aromatic_res_interface']))

############################################################################################################################
############################################################################################################################
############################################################################################################################
    with tab5:
        st.subheader("Binder Analysis Overall")

        # ensure df is defined
        df = st.session_state.df_out

        # -------------------------------
        # BAR PLOT (4Å contacts)
        # -------------------------------


        if df is not None and df.shape[0] > 0:
            fig_contacts_h = px.bar(
                df,
                x='design_id',
                y='binder_hydrophobic_res_fraction',
                title='Binder Hydrophobicity for all designs (%)'
            )
            st.plotly_chart(fig_contacts_h, use_container_width=False)
        else:
            st.info("No data to plot. Please run analysis first.")

        st.write("------------------------------------------------------------------------------------------")
        if df is not None and df.shape[0] > 0:
            fig_contacts_bc = px.bar(
                df,
                x='design_id',
                y='binder_positive_res_fraction',
                title='Binder Basic Residues for all designs (%)'
            )
            st.plotly_chart(fig_contacts_bc, use_container_width=True)
        else:
            st.info("No data to plot. Please run analysis first.")

        st.write("------------------------------------------------------------------------------------------")

        if df is not None and df.shape[0] > 0:
            fig_contacts_ac = px.bar(
                df,
                x='design_id',
                y='binder_negative_res_fraction',
                title='Binder Acidic Residues for all designs (%)'
            )
            st.plotly_chart(fig_contacts_ac, use_container_width=True)
        else:
            st.info("No data to plot. Please run analysis first.")

        st.write("------------------------------------------------------------------------------------------")


        if df is not None and df.shape[0] > 0:
            fig_contacts_ar = px.bar(
                df,
                x='design_id',
                y='binder_aromatic_res_fraction',
                title='Binder Aromatic Residues for all designs (%)'
            )
            st.plotly_chart(fig_contacts_ar, use_container_width=True)
        else:
            st.info("No data to plot. Please run analysis first.")




############################################################################################################################
############################################################################################################################
############################################################################################################################
    with tab6:

        if not check_box_family:
            st.subheader("Comming Soon!")
            pass
        else:
            st.subheader("Overall Analysis of Designs")




            def analysis_final(df_col_name):
                analysis_final = []
                for _, row in st.session_state.df_out.iterrows():
                    analysis_final.append(row[df_col_name])
                analysis_final_one_dim = []
                for i in range(len(analysis_final)):
                    for j in range(len(analysis_final[i])):
                        analysis_final_one_dim.append(analysis_final[i][j])

                n = 1
                dict_res = {}

                for i in range(len(st.session_state.df_out['target_seq'][0])):
                    dict_res[i + 1] = 0

                for i in range(len(dict_res)):
                    dict_res[i + 1] = analysis_final_one_dim.count(list(dict_res.keys())[i])

                for i in range(len(dict_res)):
                    if dict_res[i + 1] == 0:
                        del dict_res[i + 1]
                return dict_res





            analysis_final_3A = analysis_final('3A_Residues')
            analysis_final_4A = analysis_final('4A_Residues')
            analysis_final_hydrophobic = analysis_final('hydrophobic_contacts_Residues')
            analysis_final_salt_bridge = analysis_final('salt_bridge_contacts_Residues')

            target_sequence = st.session_state.df_out['target_seq'][0]
            one_to_three = {
                "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP",
                "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY",
                "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
                "M": "MET", "F": "PHE", "P": "PRO", "S": "SER",
                "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
            }
            #####################################################################################
            # Frequency of target residues at the 3 Å interface

            # Assumes target PDB numbering starts at 1.
            # add_target_res_offset makes the numbering match the original CapRel sequence.
            resname_by_resid = {
                position + 1 : one_to_three.get(aa, "UNK")
                for position, aa in enumerate(target_sequence)
            }

            # Sort residues numerically
            sorted_items_3 = sorted(
                ((int(resid), count) for resid, count in analysis_final_3A.items()),
                key=lambda x: x[0]
            )

            resids = [(resid) for resid, _ in sorted_items_3]
            counts = [count for _, count in sorted_items_3]

            residue_labels = [f"{resid}    ({resname_by_resid.get(resid, 'UNK')})"
                for resid in resids]

            fig_contacts = px.bar(
                x=residue_labels,
                y=counts,
                labels={"x": "Target residue number", "y": "Number of designs"},
                title="Frequency of target residues at the 3 Å interface"
            )

            fig_contacts.update_xaxes(
                type="category",
                tickmode="array",
                tickvals=residue_labels,
                ticktext=residue_labels,
                tickangle=-90,
            )

            fig_contacts.update_layout(
                height=600,
                margin=dict(b=140),
            )

            st.plotly_chart(fig_contacts, use_container_width=True)
             #########################################################################
            #Frequency of target residues at the 4 Å interface
            # Sort residues numerically
            sorted_items_4 = sorted(
                ((int(resid), count) for resid, count in analysis_final_4A.items()),
                key=lambda x: x[0]
            )

            resids4 = [(resid) for resid, _ in sorted_items_3]
            counts4 = [count for _, count in sorted_items_3]

            residue_labels4 = [f"{resid}    ({resname_by_resid.get(resid, 'UNK')})"
                              for resid in resids4]

            fig_contacts4 = px.bar(
                x=residue_labels4,
                y=counts4,
                labels={"x": "Target residue number", "y": "Number of designs"},
                title="Frequency of target residues at the 4 Å interface"
            )

            fig_contacts4.update_xaxes(
                type="category",
                tickmode="array",
                tickvals=residue_labels4,
                ticktext=residue_labels4,
                tickangle=-90,
            )

            fig_contacts4.update_layout(
                height=600,
                margin=dict(b=140),
            )

            st.plotly_chart(fig_contacts4, use_container_width=True)

            #st.write('Residues in 4A interface')
            #st.write(analysis_final_4A)

            #####################################################################################
            #####################################################################################
            # Frequency of target residues at the Hydrophobic contacts  4Å interface

            # Assumes target PDB numbering starts at 1.
            # add_target_res_offset makes the numbering match the original CapRel sequence.
            resname_by_resid = {
                position + 1: one_to_three.get(aa, "UNK")
                for position, aa in enumerate(target_sequence)
            }

            # Sort residues numerically
            sorted_items_hypho = sorted(
                ((int(resid), count) for resid, count in analysis_final_hydrophobic.items()),
                key=lambda x: x[0]
            )

            resids_hypho = [(resid) for resid, _ in sorted_items_hypho]
            counts_hypho = [count for _, count in sorted_items_hypho]

            residue_hypho_labels = [f"{resid}    ({resname_by_resid.get(resid, 'UNK')})"
                              for resid in resids_hypho]

            fig_contacts_hypho = px.bar(
                x=residue_hypho_labels,
                y=counts_hypho,
                labels={"x": "Target residue number", "y": "Number of designs"},
                title="Frequency of target residues at the Hydrophobic contacts at the interface"
            )

            fig_contacts_hypho.update_xaxes(
                type="category",
                tickmode="array",
                tickvals=residue_hypho_labels,
                ticktext=residue_hypho_labels,
                tickangle=-90,
            )

            fig_contacts_hypho.update_layout(
                height=600,
                margin=dict(b=140),
            )

            st.plotly_chart(fig_contacts_hypho, use_container_width=True)

            #####################################################################################
            #####################################################################################
            # Frequency of target residues at the salt bridge  4Å interface

            # Assumes target PDB numbering starts at 1.
            # add_target_res_offset makes the numbering match the original CapRel sequence.
            resname_by_resid = {
                position + 1: one_to_three.get(aa, "UNK")
                for position, aa in enumerate(target_sequence)
            }

            # Sort residues numerically
            sorted_items_salt_bridge = sorted(
                ((int(resid), count) for resid, count in analysis_final_salt_bridge.items()),
                key=lambda x: x[0]
            )

            resids_salt_bridge = [(resid) for resid, _ in sorted_items_salt_bridge]
            counts_salt_bridge = [count for _, count in sorted_items_salt_bridge]

            residue_salt_bridge_labels = [f"{resid}    ({resname_by_resid.get(resid, 'UNK')})"
                                    for resid in resids_salt_bridge]

            fig_contacts_salt_bridge = px.bar(
                x=residue_salt_bridge_labels,
                y=counts_salt_bridge,
                labels={"x": "Target residue number", "y": "Number of designs"},
                title="Frequency of target residues at the salt bridge at the interface"
            )

            fig_contacts_salt_bridge.update_xaxes(
                type="category",
                tickmode="array",
                tickvals=residue_salt_bridge_labels,
                ticktext=residue_salt_bridge_labels,
                tickangle=-90,
            )

            fig_contacts_salt_bridge.update_layout(
                height=600,
                margin=dict(b=140),
            )

            st.plotly_chart(fig_contacts_salt_bridge, use_container_width=True)

############################################################################################################################
############################################################################################################################
############################################################################################################################

    with tab7:
        st.subheader("Filter your designs")
        col1, col2 = st.columns([1, 2])
        with col1:

            st.markdown("##### Designs Filtered by Selected Metric(s)")
            filters = ['Average_pLDDT', 'n_contacts_3A_Residues', 'n_contacts_4A_Residues',
                       'Number of All H bonds','n_hydrophobic_contacts_ATOM', 'n_hydrophobic_contacts_Residues',
                       'n_salt_bridge_ATOM', 'n_salt_bridge_Residues', 'n_pi_pi_contacts', 'n_cation_pi_contacts', 'binder_hydrophobic_res_fraction',
                       'binder_hydrophobic_res_fraction_interface', 'binder_positive_res_fraction', 'binder_positive_res_fraction_interface',
                       'binder_negative_res_fraction','binder_negative_res_fraction_interface','binder_aromatic_res_fraction',
                       'binder_aromatic_res_fraction_interface']
            # cutoffs = [0.9, 0.8, 0.7, 0.5, 0.7, 0.2, 0.13, -70]
            # aa = []

            # final_filters = dict(zip(filters,cutoffs))
            df_final = st.session_state.df_out.copy()

            for filter in filters:
                show_filter = st.checkbox(filter, value=False)
                df_final[filter] = pd.to_numeric(df_final[filter], errors='coerce')
                if show_filter:
                    if ('pLDDT' in filter):
                        cutoff = st.slider("choose your cutoff (%) : ", 0, 100, 70, key=f"slider_{filter}")
                        cutoff = cutoff / 100
                        df_final = df_final[df_final[filter] >= cutoff]

                    elif ('fraction' in filter):
                        cutoff = st.slider("Choose your cutoff (%)  : ", 0, 100, 5, key=f"slider_{filter}")
                        #cutoff = cutoff / 100
                        # aa.append(cutoff)
                        df_final = df_final[df_final[filter] >= cutoff]


                    elif ('3A' in filter) or ('4A' in filter):
                        cutoff = st.slider("Choose your cutoff   : ", 0, 100, 5, key=f"slider_{filter}")
                        df_final = df_final[df_final[filter] >= cutoff]

                    elif ('pi_contacts' in filter):
                        cutoff = st.slider("Choose your cutoff   : ", 0, 10, 0, key=f"slider_{filter}")
                        df_final = df_final[df_final[filter] >= cutoff]

                    elif ('binder' in filter):
                        cutoff = st.slider("Choose your cutoff (%)   : ", 0, 100, 5, key=f"slider_{filter}")
                        df_final = df_final[df_final[filter] >= cutoff]


                    elif (filter == 'n_hydrophobic_contacts_ATOM') or (filter == 'n_hydrophobic_contacts_Residues') or (filter == 'n_salt_bridge_ATOM') or (filter ==  'n_salt_bridge_Residues'):
                        cutoff = st.slider("Choose your cutoff   : ", 0, 100, 5, key=f"slider_{filter}")
                        df_final = df_final[df_final[filter] >= cutoff]

                    elif "H bonds" in filter:
                        cutoff = st.slider("Choose your cutoff   : ", 0, 20, 1, key=f"slider_{filter}")
                        df_final = df_final[df_final[filter] >= cutoff]
                   

        with col2:
            st.markdown("##### Recommended designs based on filters: ")
            #for c in ['pairs_3A', 'pairs_4A', 'hydrophobic_contacts']:
            for c in df_final.columns:
                if c in df_final.columns:
                    df_final[c] = df_final[c].astype(str)
            # Save to session state so Streamlit tracks it
            st.session_state.df_final = df_final
            st.dataframe(df_final)
            #df_final
            if st.session_state.df_final is not None:
                st.download_button(
                    label="Download Filtered CSV",
                    data=st.session_state.df_final.to_csv(index=False),
                    file_name="Filtered_Designs.csv",
                    mime="text/csv",
                )
            #st.download_button("Download Ranked CSV", st.session_state.df_rank.to_csv(index=False), "ranked.csv")







    if st.sidebar.button(" RESET "):
        st.session_state.df_out = None
        st.session_state.df_rank = None
        st.session_state.pdb_map = {}
        st.session_state.analysis_done = False
        st.rerun()


else:
    st.info("Please upload PDB files. ")
    blankk = st.header("")
    st.markdown("#### Before starting, it is recommended to download and read the tutorial")


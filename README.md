
# PDB Analyzer Web App 

This is a Streamlit-based web application for analyzing **PDB files**, exploring **interface metrics**, and visualizing **3D protein structures**.

This tool helps you interpret the results generated after running Protein Design Software like AlphaFold, BindCraft, including the Accepted PDB files.




---

## Getting Started
### Accessing the PDB Analyzer Web App
This application is hosted on Streamlit Cloud, providing easy access without any local installation required. The link :
https://pdbanalyzermjshadfar.streamlit.app/

 
(Includes a RESET button in the sidebar to clear all information.)

## 🚀 Features

### ✅ 1. Summary Tab
- Displays all analyzed results in a single interactive table  
  
- Provides buttons to download:
  - summary CSV  
  - Ranked designs CSV  

---

### ✅ 2. Visualization Tab
- Interactive **3D structure viewer** using *py3Dmol*
- Customizable chain coloring  
- Counts and plots the number of **3 Å and 4 Å contacts** between binder and target  

---

### ✅ 3. Details Tab
For each design, you can inspect:
- Target and binder sequences  
- Residue–residue contact pairs (3 Å and 4 Å)   


---

## 📁 Input Requirements

Upload the following:

1. **Accepted PDB designs**  
   Output from BindCraft or AlphaFold.



The app will automatically match each PDB to its corresponding row in the CSV.

---

## 🛠️ How to use it?

You have two choices:  
- install all the dependencies.
- just use the Wep App via the link.  (https://bindcraftoutputanalyse-ednrdmyx97eon4hlfkjatb.streamlit.app/)


If you want to install it locally:
####  Dependencies 

git clone https://github.com/mj72git/PDB-Analyzer.git
cd PDB-Analyzer-main

- Python
- MDAnalysis
- numpy
- pandas
- plotly
- py3Dmol
- streamlit





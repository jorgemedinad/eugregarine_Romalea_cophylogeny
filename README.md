# README

### **Authors**  
Jorge H. Medina-Duran

---

### **Project Title**  
**Colonization and Vicariance Structure Eugregarine Associations in *Romalea* Grasshoppers**

---

### **Description**  
This repository contains all scripts, datasets, analyses, and outputs required to reproduce the event-based and global-fit cophylogenetic analyses presented in the study *"Cophylogenetic Analysis of Eugregarine Parasites in Romalea Grasshoppers"*.  

The project integrates event-based reconciliation (eMPRess) and global-fit methods (PACo, ParaFit) to assess cospeciation and host-switching dynamics between *Romalea* hosts and their eugregarine parasites.

---

### **Directory Structure**
```
.
├── input/                         # Input files: occurrence matrices, phylogenies, collecting data, shapefiles
├── output/                        # All generated outputs: raw figures, Empress results, PACo/ParaFit analyses
├── scripts/                       # All R and bash scripts organized by analysis step
├── manuscript_plates/             # Paper figures
├── supplementary_information/     # All supplementary figures and tables for the manuscript
├── README.md                      # Project documentation
```

---

### **Pipeline Overview and Script Descriptions**

| Script Name                   | Description |
|:-------------------------------|:------------|
| **01_prepare_occurrence_matrix.R** | Deduplicates FASTA sequences, matches host collection data, builds parasite-host occurrence matrix. |
| **02_trim_trees_occurrence_matrix.R** | Trims host and parasite trees and the occurrence matrix to ensure matching taxa before analyses. |
| **03_subsample_one_to_one_matrices.R** | Generates 50 one-to-one parasite-host interaction matrices for event-based reconciliation. |
| **04_run_empress.sh**         | Automates Empress reconciliation and p-value calculations across cost schemes and replicates. |
| **05_analyze_empress_results.R** | Summarizes Empress outputs: event counts, p-values, and visualizations of cost scheme results. |
| **06_paco_parafit_analysis.R** | Performs PACo and ParaFit global-fit analyses, computes link contributions and significance. |
| **07_map_collecting_sites.R** | Maps Romalea collecting sites, overlays them on elevation and biogeographic province maps. |

---

### **Running Empress Analyses**

The `04_run_empress.sh` script runs Empress reconciliations automatically across cost schemes.  

**Step 1: Create and activate Python environment**
```bash
conda create -n empress_env python=3.7
conda activate empress_env
```
Install eMPRess following the [official instructions](https://github.com/ssantichaivekin/empress/wiki/Install-Empress-with-Command-Line-Interface-for-Development).  
You may also need:

```bash
conda install matplotlib biopython shapely networkx
```

**Step 2: Configure script paths**  
Edit `04_run_empress.sh`:
- `PYTHON_EXEC` → path to the Python binary inside your conda environment
- `data_path` → path to project `output/` folder

**Step 3: Run the script**  
Place the script inside your Empress installation (`empress/` folder) and execute:
```bash
bash 04_run_empress.sh
```

---

### **How to Reproduce the Project**

1. **Update working paths**  
   Edit the top of each script to define `input/` and `output/` directories.

2. **Install required R packages**
```r
install.packages(c("ape", "phytools", "paco", "ggplot2", "tidyr", "sf", 
                   "rnaturalearth", "raster", "geodata", "terra", "viridis"))
```

3. **Run analysis pipeline sequentially**  
   - Execute scripts `01_prepare_occurrence_matrix.R` → `07_map_collecting_sites.R` in order.
   - Execute `04_run_empress.sh` separately under Linux/WSL for Empress reconciliations.

4. **Output locations**  
   - Raw figures saved in `output/figures/`
   - Empress results in `output/results_empress/`
   - PACo/ParaFit summaries in `output/paco/` and `output/parafit/`
   - Province assignments in `output/geographic_provinces/`

---
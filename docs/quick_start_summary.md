# Quick Start: Pre-equilibrated Membrane + PPM 3.0 Workflow

## 🚀 **TL;DR - Essential Steps**

### **1. Extract Membrane Frame**
```bash
echo 0 | gmx trjconv -f membrane_equilibration.xtc -s membrane.tpr -o membrane_final.gro -dump 1000
```

### **2. Position Protein with PPM 3.0**
- **Web**: Go to https://opm.phar.umich.edu/ppm_server
- **Upload**: Your protein PDB file
- **Download**: Positioned protein coordinates

### **3. Combine Systems**
```python
import MDAnalysis as mda
protein = mda.Universe("positioned_protein.pdb")
membrane = mda.Universe("membrane_final.gro") 
combined = mda.Merge(protein.atoms, membrane.atoms)
combined.atoms.write("combined_system.gro")
```

### **4. Setup GROMACS Simulation**
```bash
# Generate protein topology
gmx pdb2gmx -f positioned_protein.pdb -o protein.gro -p protein.top

# Combine topologies (manual editing needed)
# Create equilibration.mdp file

# Run simulation
gmx grompp -f equilibration.mdp -c combined_system.gro -p combined_system.top -o eq.tpr
gmx mdrun -deffnm equilibration
```

### **5. Analyze with OrientationAnalyzer**
```python
from memprot.analyzers import OrientationAnalyzer
# Configure, run analysis, get Dz and θ parameters
```

---

## 📁 **Required Files**

### **Input:**
- `membrane_equilibration.xtc` - Your GROMACS membrane trajectory
- `membrane.tpr` - GROMACS run input file
- `protein.pdb` - Your protein structure

### **Generated:**
- `membrane_final.gro` - Extracted membrane frame
- `positioned_protein.pdb` - PPM 3.0 positioned protein
- `combined_system.gro` - Final system for simulation
- `combined_system.top` - Topology file

---

## 🔧 **Key Tools**

1. **PPM 3.0**: Protein positioning in membranes
2. **GROMACS**: MD simulation engine  
3. **MDAnalysis**: System manipulation
4. **OrientationAnalyzer**: K-Ras4B style analysis

---

## ⚠️ **Common Pitfalls**

- **Topology mismatch**: Ensure protein and membrane force fields are compatible
- **Box dimensions**: Protein must fit within membrane box
- **Equilibration**: Use position restraints initially
- **PPM 3.0 timeout**: Have manual positioning as backup

---

## 📊 **Analysis Output**

The `OrientationAnalyzer` provides:
- **Dz**: Distance between effector lobe COM and membrane COM (Z-direction)
- **θ (theta)**: Tilt angle of protein relative to membrane normal
- **Time series plots**: Evolution of orientation parameters
- **Statistical summaries**: Mean, std, min, max values

---

## 🔗 **Resources**

- **PPM 3.0**: https://opm.phar.umich.edu/ppm_server
- **CHARMM-GUI**: http://www.charmm-gui.org/ (alternative approach)
- **Full Guide**: `membrane_protein_setup_guide.md`
- **Workflow Script**: `membrane_protein_setup_workflow.py`

This workflow bridges the gap between your pre-equilibrated GROMACS membrane and the powerful orientation analysis capabilities of the `memprot` package! 
# Guide: Using Pre-equilibrated GROMACS Membrane with PPM 3.0 for Protein Positioning

This guide walks you through the process of taking a pre-equilibrated membrane from GROMACS and positioning a protein on its surface using the PPM 3.0 service.

## Overview

The workflow involves:
1. **Membrane Preparation**: Extract and analyze your pre-equilibrated GROMACS membrane
2. **Protein Positioning**: Use PPM 3.0 to position your protein in/on the membrane
3. **System Combination**: Merge the positioned protein with your membrane
4. **Simulation Setup**: Prepare files for MD simulation
5. **Analysis**: Use the `OrientationAnalyzer` to study the results

---

## Step 1: Prepare Your Pre-equilibrated Membrane

### 1.1 Extract Membrane Frame
From your GROMACS membrane equilibration trajectory, extract a representative frame:

```bash
# Extract the last frame of equilibration
echo 0 | gmx trjconv -f membrane_equilibration.xtc -s membrane.tpr -o membrane_final.gro -dump 1000

# Or extract a specific frame
echo 0 | gmx trjconv -f membrane_equilibration.xtc -s membrane.tpr -o membrane_final.gro -b 1000 -e 1000
```

### 1.2 Analyze Membrane Properties
Use MDAnalysis to understand your membrane:

```python
import MDAnalysis as mda
import numpy as np

# Load membrane
u = mda.Universe("membrane_final.gro")

# Find phosphate groups (adjust selection for your lipids)
phosphates = u.select_atoms("name P or name PO4")

# Get membrane properties
z_coords = phosphates.positions[:, 2]
membrane_center = np.mean(z_coords)
membrane_thickness = np.max(z_coords) - np.min(z_coords)

print(f"Membrane center Z: {membrane_center:.2f} Å")
print(f"Membrane thickness: {membrane_thickness:.2f} Å")
print(f"Box dimensions: {u.dimensions}")
```

---

## Step 2: Position Protein Using PPM 3.0

### 2.1 PPM 3.0 Web Server (Recommended)

1. **Go to the PPM server**: https://opm.phar.umich.edu/ppm_server
2. **Upload your protein PDB file** or paste the content
3. **Select calculation type**: 
   - "Position" for membrane positioning
   - "Energy" if you only want binding energy
4. **Choose membrane type**: 
   - DOPC (default)
   - Or select the lipid type matching your membrane
5. **Submit and download results**

### 2.2 Alternative: Manual PPM 3.0 API Call

```python
import requests

def use_ppm_server(protein_pdb_path, output_path="positioned_protein.pdb"):
    """Submit protein to PPM 3.0 server for positioning."""
    
    # Read protein PDB
    with open(protein_pdb_path, 'r') as f:
        pdb_content = f.read()
    
    # PPM server endpoint
    url = "https://opm.phar.umich.edu/ppm_server"
    
    # Prepare data
    data = {
        'pdb_content': pdb_content,
        'membrane_type': 'DOPC',  # Adjust as needed
        'calculation_type': 'position'
    }
    
    # Submit request
    response = requests.post(url, data=data, timeout=300)
    
    if response.status_code == 200:
        # Parse and save positioned protein
        positioned_pdb = parse_ppm_response(response.text)
        with open(output_path, 'w') as f:
            f.write(positioned_pdb)
        print(f"Positioned protein saved to {output_path}")
        return output_path
    else:
        raise Exception(f"PPM server error: {response.status_code}")

# Usage
positioned_protein = use_ppm_server("your_protein.pdb")
```

### 2.3 Fallback: Manual Positioning

If PPM 3.0 is unavailable, you can manually position the protein:

```python
import MDAnalysis as mda
import numpy as np

def manual_position_protein(protein_pdb, membrane_gro, output_pdb, surface_offset=10.0):
    """Manually position protein at membrane surface."""
    
    # Load systems
    protein = mda.Universe(protein_pdb)
    membrane = mda.Universe(membrane_gro)
    
    # Get membrane surface Z-coordinate
    phosphates = membrane.select_atoms("name P or name PO4")
    upper_surface_z = np.max(phosphates.positions[:, 2])
    
    # Position protein above membrane surface
    protein_com = protein.atoms.center_of_mass()
    target_z = upper_surface_z + surface_offset
    
    # Translate protein
    translation = np.array([0, 0, target_z - protein_com[2]])
    protein.atoms.translate(translation)
    
    # Center in XY plane
    membrane_center_xy = membrane.dimensions[:2] / 2
    protein_center_xy = protein.atoms.center_of_mass()[:2]
    xy_translation = np.array([membrane_center_xy[0] - protein_center_xy[0],
                              membrane_center_xy[1] - protein_center_xy[1], 0])
    protein.atoms.translate(xy_translation)
    
    # Save positioned protein
    protein.atoms.write(output_pdb)
    print(f"Manually positioned protein saved to {output_pdb}")

# Usage
manual_position_protein("protein.pdb", "membrane_final.gro", "positioned_protein.pdb")
```

---

## Step 3: Combine Protein with Membrane

### 3.1 Using MDAnalysis

```python
import MDAnalysis as mda

def combine_systems(positioned_protein_pdb, membrane_gro, output_gro="combined_system.gro"):
    """Combine positioned protein with membrane."""
    
    # Load systems
    protein = mda.Universe(positioned_protein_pdb)
    membrane = mda.Universe(membrane_gro)
    
    # Merge systems
    combined = mda.Merge(protein.atoms, membrane.atoms)
    combined.dimensions = membrane.dimensions  # Use membrane box dimensions
    
    # Write combined system
    combined.atoms.write(output_gro)
    print(f"Combined system written to {output_gro}")
    
    return output_gro

# Usage
combined_system = combine_systems("positioned_protein.pdb", "membrane_final.gro")
```

### 3.2 Using GROMACS Tools

```bash
# Convert protein PDB to GRO with correct box
gmx editconf -f positioned_protein.pdb -o protein.gro -box $(gmx dump -s membrane.tpr | grep box | tail -1)

# Combine protein and membrane
gmx insert-molecules -f membrane_final.gro -ci protein.gro -o combined_system.gro -nmol 1
```

---

## Step 4: Setup Topology and Simulation Files

### 4.1 Create Combined Topology

```bash
# Create combined topology file
cat > combined_system.top << 'EOF'
; Combined protein-membrane topology
#include "amber99sb-ildn.ff/forcefield.itp"

; Include protein topology (generated with pdb2gmx)
#include "protein.itp"

; Include membrane lipid topologies
#include "lipids.itp"

; Include water model
#include "amber99sb-ildn.ff/tip3p.itp"

[ system ]
Protein in membrane

[ molecules ]
; Molecule name    Number
Protein            1
POPC              128    ; Adjust numbers based on your membrane
POPE               32
SOL               5000   ; Adjust water molecules
EOF
```

### 4.2 Generate Protein Topology

```bash
# Generate protein topology with pdb2gmx
gmx pdb2gmx -f positioned_protein.pdb -o protein_processed.gro -p protein.top -ignh

# Extract protein.itp from the generated topology
```

### 4.3 Create MDP File for Equilibration

```bash
cat > equilibration.mdp << 'EOF'
; Membrane-protein equilibration parameters
integrator = md
dt = 0.002
nsteps = 250000  ; 500 ps

; Output control
nstxout = 5000
nstvout = 5000
nstenergy = 1000
nstlog = 1000

; Temperature coupling
tcoupl = V-rescale
tc-grps = Protein Membrane SOL
tau-t = 0.1 0.1 0.1
ref-t = 310 310 310

; Pressure coupling
pcoupl = Parrinello-Rahman
pcoupltype = semiisotropic
tau-p = 2.0
ref-p = 1.0 1.0
compressibility = 4.5e-5 4.5e-5

; Constraints
constraints = h-bonds
constraint-algorithm = LINCS

; Electrostatics
coulombtype = PME
rcoulomb = 1.2

; Van der Waals
vdwtype = Cut-off
rvdw = 1.2

; PBC
pbc = xyz

; Position restraints (optional for initial equilibration)
define = -DPOSRES
EOF
```

---

## Step 5: Run Simulation

### 5.1 Preprocessing

```bash
# Generate run input file
gmx grompp -f equilibration.mdp -c combined_system.gro -p combined_system.top -o equilibration.tpr

# Check for warnings and fix topology if needed
```

### 5.2 Run Equilibration

```bash
# Run equilibration
gmx mdrun -deffnm equilibration -v

# Monitor progress
tail -f equilibration.log
```

### 5.3 Production Run

```bash
# Create production MDP (longer simulation, less frequent output)
# Then run production simulation
gmx grompp -f production.mdp -c equilibration.gro -p combined_system.top -o production.tpr
gmx mdrun -deffnm production -v
```

---

## Step 6: Analysis with OrientationAnalyzer

Once your simulation is complete, use the `OrientationAnalyzer` to study protein-membrane interactions:

```python
from memprot.core.config import AnalysisConfig
from memprot.core.system import SimulationSystem
from memprot.analyzers import OrientationAnalyzer

# Setup analysis configuration
config = AnalysisConfig(
    start_time=100000,  # Skip first 100 ns for equilibration
    dt=10,
    output_dir="orientation_analysis"
)

# Create simulation system
system = SimulationSystem(
    name="protein_membrane",
    topology_file="combined_system.gro",
    trajectory_file="production.xtc",
    protein_selection="protein",
    membrane_selection="resname POPC POPE",  # Adjust for your lipids
)

# Initialize orientation analyzer
analyzer = OrientationAnalyzer(
    config=config,
    effector_lobe_residues="resid 1-50",     # Adjust for your protein
    allosteric_lobe_residues="resid 51-100", # Adjust for your protein
    beta4_residues="resid 30-35",            # Key secondary structure
    beta5_residues="resid 60-65"             # Key secondary structure
)

# Add system and run analysis
analyzer.add_system(system)
analyzer.run_analysis()

# Get results
summary = analyzer.get_orientation_summary("protein_membrane")
print(f"Mean Dz: {summary['mean_dz']:.2f} ± {summary['std_dz']:.2f} Å")
print(f"Mean θ: {summary['mean_theta']:.2f} ± {summary['std_theta']:.2f}°")

# Plot time series
analyzer.plot_orientation_timeseries("protein_membrane", "orientation_analysis.png")

# Save results
analyzer.save_results()
```

---

## Tips and Troubleshooting

### Common Issues

1. **PPM 3.0 Server Unavailable**
   - Use the manual positioning approach
   - Check if protein structure is valid PDB format
   - Try smaller protein structures first

2. **Topology Errors**
   - Ensure protein force field matches membrane force field
   - Check that all residues are recognized by `pdb2gmx`
   - Add missing lipid parameters to topology

3. **System Crashes During Simulation**
   - Use position restraints during initial equilibration
   - Check for clashes between protein and membrane
   - Start with shorter time steps (dt = 0.001)

### Best Practices

1. **Equilibration Strategy**
   - Start with position restraints on protein
   - Gradually remove restraints over multiple equilibration stages
   - Monitor system energy and temperature

2. **Membrane Selection**
   - Use realistic lipid compositions for your system
   - Consider asymmetric membranes if biologically relevant
   - Ensure adequate membrane size to avoid periodic interactions

3. **Analysis Considerations**
   - Allow sufficient equilibration time (>100 ns)
   - Verify protein doesn't unfold during simulation
   - Check membrane integrity throughout simulation

---

## Alternative Approaches

### Using CHARMM-GUI

If you prefer a web-based approach, you can also use CHARMM-GUI:

1. Go to http://www.charmm-gui.org/
2. Use "Membrane Builder" → "Bilayer Builder"
3. Upload your positioned protein
4. Select membrane composition matching your pre-equilibrated system
5. Follow the CHARMM-GUI workflow

### Using Other Positioning Tools

- **MEMPROT**: For membrane protein orientation prediction
- **MemProtMD**: Database of pre-simulated membrane proteins
- **Rosetta MP**: For membrane protein modeling and design

---

This workflow provides a comprehensive approach to combining pre-equilibrated membranes with protein positioning using PPM 3.0, setting up simulations, and analyzing the results with the `OrientationAnalyzer` tool. 
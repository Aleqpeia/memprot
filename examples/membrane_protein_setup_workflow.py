#!/usr/bin/env python3
"""
Workflow for setting up protein-membrane systems using pre-equilibrated membranes and PPM 3.0.

This script demonstrates how to:
1. Process pre-equilibrated GROMACS membrane
2. Use PPM 3.0 for protein positioning
3. Combine systems for MD simulation
4. Analyze with the memprot OrientationAnalyzer

Requirements:
- MDAnalysis
- requests (for PPM 3.0 API)
- memprot package (for analysis)
"""

import MDAnalysis as mda
import numpy as np
import requests
import subprocess
import tempfile
import os
from pathlib import Path
from typing import Optional, Tuple, Dict
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class MembraneProteinSetup:
    """Class to handle membrane-protein system setup."""
    
    def __init__(self, gromacs_membrane_gro: str, gromacs_membrane_top: str):
        """Initialize with pre-equilibrated GROMACS membrane.
        
        Args:
            gromacs_membrane_gro: Path to equilibrated membrane .gro file
            gromacs_membrane_top: Path to membrane .top file
        """
        self.membrane_gro = gromacs_membrane_gro
        self.membrane_top = gromacs_membrane_top
        self.membrane_universe = None
        self.protein_positioned = None
        
    def load_membrane(self) -> None:
        """Load the pre-equilibrated membrane system."""
        logger.info("Loading pre-equilibrated membrane from GROMACS...")
        try:
            self.membrane_universe = mda.Universe(self.membrane_gro)
            logger.info(f"Loaded membrane with {len(self.membrane_universe.atoms)} atoms")
            logger.info(f"Box dimensions: {self.membrane_universe.dimensions}")
            
            # Get membrane center and thickness
            self._analyze_membrane_properties()
            
        except Exception as e:
            logger.error(f"Error loading membrane: {e}")
            raise
    
    def _analyze_membrane_properties(self) -> None:
        """Analyze membrane properties for positioning."""
        # Find lipid phosphate groups (common reference)
        try:
            phosphates = self.membrane_universe.select_atoms("name P or name PO4")
            if len(phosphates) == 0:
                # Fallback: try other lipid head group atoms
                phosphates = self.membrane_universe.select_atoms("name N")
                
            if len(phosphates) > 0:
                z_coords = phosphates.positions[:, 2]
                self.membrane_center_z = np.mean(z_coords)
                self.membrane_thickness = np.max(z_coords) - np.min(z_coords)
                
                # Identify upper and lower leaflets
                upper_leaflet = phosphates[z_coords > self.membrane_center_z]
                lower_leaflet = phosphates[z_coords <= self.membrane_center_z]
                
                self.upper_leaflet_z = np.mean(upper_leaflet.positions[:, 2])
                self.lower_leaflet_z = np.mean(lower_leaflet.positions[:, 2])
                
                logger.info(f"Membrane center Z: {self.membrane_center_z:.2f} Å")
                logger.info(f"Membrane thickness: {self.membrane_thickness:.2f} Å")
                logger.info(f"Upper leaflet Z: {self.upper_leaflet_z:.2f} Å")
                logger.info(f"Lower leaflet Z: {self.lower_leaflet_z:.2f} Å")
            else:
                logger.warning("Could not find phosphate groups for membrane analysis")
                self.membrane_center_z = self.membrane_universe.dimensions[2] / 2
                
        except Exception as e:
            logger.warning(f"Membrane analysis failed: {e}")
            self.membrane_center_z = self.membrane_universe.dimensions[2] / 2

    def use_ppm_server(self, protein_pdb: str, output_file: str = "protein_positioned.pdb") -> str:
        """Use PPM 3.0 web server to position protein in membrane.
        
        Args:
            protein_pdb: Path to protein PDB file
            output_file: Output file for positioned protein
            
        Returns:
            Path to positioned protein file
        """
        logger.info("Using PPM 3.0 server for protein positioning...")
        
        # PPM 3.0 server URL
        ppm_url = "https://opm.phar.umich.edu/ppm_server"
        
        try:
            # Read protein PDB file
            with open(protein_pdb, 'r') as f:
                pdb_content = f.read()
            
            # Prepare request data
            data = {
                'pdb_text': pdb_content,
                'calculation_type': 'positioning',  # or 'energy' for binding energy only
                'membrane_type': 'DOPC',  # Default membrane type
            }
            
            logger.info("Submitting to PPM 3.0 server...")
            response = requests.post(ppm_url, data=data, timeout=300)
            
            if response.status_code == 200:
                # Parse response to get positioned protein
                # Note: Actual PPM server response format may vary
                # This is a simplified example
                positioned_pdb = self._parse_ppm_response(response.text)
                
                with open(output_file, 'w') as f:
                    f.write(positioned_pdb)
                    
                logger.info(f"Positioned protein saved to {output_file}")
                return output_file
            else:
                logger.error(f"PPM server request failed: {response.status_code}")
                raise Exception(f"PPM server error: {response.status_code}")
                
        except Exception as e:
            logger.error(f"Error using PPM server: {e}")
            # Fallback: manual positioning
            return self._manual_protein_positioning(protein_pdb, output_file)
    
    def _parse_ppm_response(self, response_text: str) -> str:
        """Parse PPM server response to extract positioned protein.
        
        Note: This is a simplified parser. The actual PPM response format
        may require more sophisticated parsing.
        """
        # PPM typically returns the positioned protein coordinates
        # along with membrane boundary information
        lines = response_text.split('\n')
        pdb_lines = []
        
        for line in lines:
            if line.startswith(('ATOM', 'HETATM', 'TER', 'END')):
                pdb_lines.append(line)
        
        return '\n'.join(pdb_lines)
    
    def _manual_protein_positioning(self, protein_pdb: str, output_file: str) -> str:
        """Fallback manual protein positioning based on membrane properties.
        
        This positions the protein at the membrane surface based on the
        pre-equilibrated membrane properties.
        """
        logger.info("Using manual positioning as fallback...")
        
        # Load protein
        protein_u = mda.Universe(protein_pdb)
        
        # Get protein center of mass
        protein_com = protein_u.atoms.center_of_mass()
        
        # Position protein at membrane surface
        # Choose upper or lower leaflet (here we choose upper)
        target_z = self.upper_leaflet_z + 10.0  # 10 Å above membrane surface
        
        # Calculate translation needed
        translation = np.array([0, 0, target_z - protein_com[2]])
        
        # Translate protein
        protein_u.atoms.translate(translation)
        
        # Center protein in XY plane of membrane box
        membrane_center_xy = [self.membrane_universe.dimensions[0]/2, 
                             self.membrane_universe.dimensions[1]/2]
        protein_center_xy = protein_u.atoms.center_of_mass()[:2]
        
        xy_translation = np.array([membrane_center_xy[0] - protein_center_xy[0],
                                  membrane_center_xy[1] - protein_center_xy[1], 0])
        
        protein_u.atoms.translate(xy_translation)
        
        # Save positioned protein
        protein_u.atoms.write(output_file)
        logger.info(f"Manually positioned protein saved to {output_file}")
        
        return output_file
    
    def combine_systems(self, positioned_protein_pdb: str, output_prefix: str = "combined_system") -> Tuple[str, str]:
        """Combine positioned protein with pre-equilibrated membrane.
        
        Args:
            positioned_protein_pdb: Path to PPM-positioned protein
            output_prefix: Prefix for output files
            
        Returns:
            Tuple of (combined_gro_file, combined_top_file)
        """
        logger.info("Combining protein with membrane...")
        
        # Load positioned protein
        protein_u = mda.Universe(positioned_protein_pdb)
        
        # Create combined system
        combined_gro = f"{output_prefix}.gro"
        combined_top = f"{output_prefix}.top"
        
        try:
            # Method 1: Use MDAnalysis to merge systems
            all_atoms = mda.Merge(protein_u.atoms, self.membrane_universe.atoms)
            
            # Update box dimensions to match membrane
            all_atoms.dimensions = self.membrane_universe.dimensions
            
            # Write combined structure
            all_atoms.atoms.write(combined_gro)
            
            # Handle topology file (simplified)
            self._create_combined_topology(combined_top, positioned_protein_pdb)
            
            logger.info(f"Combined system written to {combined_gro}")
            return combined_gro, combined_top
            
        except Exception as e:
            logger.error(f"Error combining systems: {e}")
            # Alternative: use GROMACS tools if available
            return self._gromacs_combine_systems(positioned_protein_pdb, output_prefix)
    
    def _create_combined_topology(self, output_top: str, protein_pdb: str) -> None:
        """Create a basic combined topology file.
        
        Note: This is simplified. For production simulations, you'll need
        proper force field parameters for the protein.
        """
        logger.info("Creating combined topology file...")
        
        with open(output_top, 'w') as f:
            f.write("; Combined protein-membrane topology\n")
            f.write("; Generated by membrane_protein_setup_workflow.py\n")
            f.write("\n")
            f.write("#include \"protein.itp\"  ; Include protein topology\n")
            f.write("#include \"membrane.itp\" ; Include membrane topology\n")
            f.write("\n")
            f.write("[ system ]\n")
            f.write("Protein in membrane\n")
            f.write("\n")
            f.write("[ molecules ]\n")
            f.write("Protein  1\n")
            f.write("; Add membrane molecule counts here\n")
        
        logger.info(f"Basic topology written to {output_top}")
        logger.warning("Please update topology with proper protein and membrane parameters!")
    
    def _gromacs_combine_systems(self, protein_pdb: str, output_prefix: str) -> Tuple[str, str]:
        """Use GROMACS tools to combine systems (if available)."""
        logger.info("Attempting to use GROMACS tools for system combination...")
        
        try:
            # Convert protein PDB to GRO
            protein_gro = f"{output_prefix}_protein.gro"
            subprocess.run(['editconf', '-f', protein_pdb, '-o', protein_gro], 
                          check=True, capture_output=True)
            
            # Use genconf or similar to combine
            combined_gro = f"{output_prefix}.gro"
            # This is a simplified example - actual GROMACS workflow may vary
            subprocess.run(['genconf', '-f', self.membrane_gro, '-o', combined_gro], 
                          check=True, capture_output=True)
            
            return combined_gro, f"{output_prefix}.top"
            
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            logger.error(f"GROMACS tools not available or failed: {e}")
            raise Exception("Both MDAnalysis and GROMACS methods failed")

    def setup_md_simulation(self, combined_gro: str, output_dir: str = "md_setup") -> Dict[str, str]:
        """Set up MD simulation files.
        
        Args:
            combined_gro: Path to combined system structure
            output_dir: Directory for MD setup files
            
        Returns:
            Dictionary with paths to setup files
        """
        logger.info("Setting up MD simulation files...")
        
        os.makedirs(output_dir, exist_ok=True)
        
        files = {}
        
        # Create basic MDP file for equilibration
        mdp_file = f"{output_dir}/equilibration.mdp"
        self._create_mdp_file(mdp_file)
        files['mdp'] = mdp_file
        
        # Copy structure file
        structure_file = f"{output_dir}/system.gro"
        subprocess.run(['cp', combined_gro, structure_file], check=True)
        files['structure'] = structure_file
        
        logger.info(f"MD setup files created in {output_dir}")
        return files
    
    def _create_mdp_file(self, mdp_file: str) -> None:
        """Create a basic MDP file for membrane-protein equilibration."""
        mdp_content = """
; Basic equilibration MDP for membrane-protein system
; Adjust parameters as needed for your specific system

integrator = md
dt = 0.002
nsteps = 50000  ; 100 ps

; Output control
nstxout = 1000
nstvout = 1000
nstenergy = 100
nstlog = 100

; Temperature coupling
tcoupl = V-rescale
tc-grps = Protein Membrane
tau-t = 0.1 0.1
ref-t = 310 310

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

; Freeze groups (optional - freeze membrane during initial equilibration)
; freezegrps = Membrane
; freezedim = Y Y Y
"""
        
        with open(mdp_file, 'w') as f:
            f.write(mdp_content)


def main():
    """Example usage of the membrane-protein setup workflow."""
    
    # Example file paths - replace with your actual files
    membrane_gro = "path/to/your/equilibrated_membrane.gro"
    membrane_top = "path/to/your/membrane.top"
    protein_pdb = "path/to/your/protein.pdb"
    
    # Check if files exist
    if not all(os.path.exists(f) for f in [membrane_gro, protein_pdb]):
        logger.error("Required input files not found. Please check paths.")
        print("\nExample usage:")
        print("1. Prepare your pre-equilibrated GROMACS membrane (.gro and .top files)")
        print("2. Have your protein structure (.pdb file)")
        print("3. Update the file paths in this script")
        print("4. Run the workflow")
        return
    
    # Initialize setup
    setup = MembraneProteinSetup(membrane_gro, membrane_top)
    
    try:
        # Step 1: Load and analyze membrane
        setup.load_membrane()
        
        # Step 2: Position protein using PPM 3.0
        positioned_protein = setup.use_ppm_server(protein_pdb, "positioned_protein.pdb")
        
        # Step 3: Combine systems
        combined_gro, combined_top = setup.combine_systems(positioned_protein, "protein_membrane_system")
        
        # Step 4: Setup MD simulation
        md_files = setup.setup_md_simulation(combined_gro, "md_simulation_setup")
        
        print(f"\n✅ Workflow completed successfully!")
        print(f"Combined system: {combined_gro}")
        print(f"MD setup directory: md_simulation_setup/")
        print(f"\nNext steps:")
        print(f"1. Review and adjust the topology file: {combined_top}")
        print(f"2. Check the MDP parameters: {md_files['mdp']}")
        print(f"3. Run GROMACS preprocessing: grompp -f {md_files['mdp']} -c {md_files['structure']} -p {combined_top}")
        print(f"4. Run simulation: mdrun")
        print(f"5. Analyze with memprot OrientationAnalyzer!")
        
    except Exception as e:
        logger.error(f"Workflow failed: {e}")
        print(f"\n❌ Workflow failed. Check the logs above.")


if __name__ == "__main__":
    main() 
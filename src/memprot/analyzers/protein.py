"""Protein structure and dynamics analysis."""

from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
from tqdm import tqdm
from loguru import logger
import MDAnalysis as mda
from MDAnalysis.analysis import rms, contacts, secondary_structure

from ..core.base import BaseAnalyzer


class ProteinAnalyzer(BaseAnalyzer):
    """Analyzes protein structure and dynamics."""
    
    def run_analysis(self, system_names: Optional[List[str]] = None) -> None:
        """Run protein analysis.
        
        Args:
            system_names: Names of systems to analyze. If None, analyze all.
        """
        if system_names is None:
            system_names = [name for name, sys in self._systems.items() if sys.has_protein]
        else:
            # Filter to only systems with protein
            system_names = [name for name in system_names 
                           if name in self._systems and self._systems[name].has_protein]
            
        if not system_names:
            logger.warning("No systems with protein found")
            return
            
        self._validate_systems(system_names)
        
        for system_name in system_names:
            logger.info(f"Analyzing protein in system: {system_name}")
            
            # Basic structural properties
            self._analyze_structure(system_name)
            
            # RMSD analysis
            self._analyze_rmsd(system_name)
            
            # Secondary structure
            self._analyze_secondary_structure(system_name)
            
            # Radius of gyration
            self._analyze_radius_of_gyration(system_name)
            
            # End-to-end distance
            self._analyze_end_to_end_distance(system_name)
            
            # Contact maps
            self._analyze_contact_maps(system_name)
            
        logger.info("Protein analysis completed")
    
    def _analyze_structure(self, system_name: str) -> None:
        """Analyze basic protein structure."""
        system = self.get_system(system_name)
        protein = system.protein
        
        if protein is None:
            return
            
        # Basic information
        n_residues = len(protein.residues)
        n_atoms = len(protein.atoms)
        
        # Residue types
        residue_types = [res.resname for res in protein.residues]
        residue_counts = pd.Series(residue_types).value_counts()
        
        # Store results
        if "structure_info" not in self.results:
            self.results["structure_info"] = {}
            
        self.results["structure_info"][system_name] = {
            "n_residues": n_residues,
            "n_atoms": n_atoms,
            "residue_counts": residue_counts.to_dict(),
            "sequence": "".join([res.resname for res in protein.residues])
        }
        
        logger.info(f"System {system_name}: {n_residues} residues, {n_atoms} atoms")
    
    def _analyze_rmsd(self, system_name: str) -> None:
        """Calculate RMSD trajectory."""
        system = self.get_system(system_name)
        protein = system.protein
        
        if protein is None:
            return
            
        # Align trajectory first
        system.align_protein_trajectory()
        
        # Calculate RMSD for backbone atoms
        ca_atoms = protein.select_atoms("name CA")
        bb_atoms = protein.select_atoms("backbone")
        
        frame_slice = self.get_frame_slice(system_name)
        
        # Set reference frame
        system.universe.trajectory[frame_slice[0]]
        ref_ca = ca_atoms.positions.copy()
        ref_bb = bb_atoms.positions.copy()
        
        rmsd_ca = []
        rmsd_bb = []
        
        for frame_idx in tqdm(frame_slice, desc=f"RMSD {system_name}"):
            system.universe.trajectory[frame_idx]
            
            # CA RMSD
            rmsd_ca.append(rms.rmsd(ca_atoms.positions, ref_ca))
            
            # Backbone RMSD
            rmsd_bb.append(rms.rmsd(bb_atoms.positions, ref_bb))
        
        # Create DataFrame
        data = {
            "rmsd_ca": np.array(rmsd_ca),
            "rmsd_backbone": np.array(rmsd_bb)
        }
        
        df = self._create_dataframe(data, system_name)
        
        if "rmsd" not in self.results:
            self.results["rmsd"] = {}
        self.results["rmsd"][system_name] = df
        
        logger.info(f"RMSD analysis completed for {system_name}")
    
    def _analyze_secondary_structure(self, system_name: str) -> None:
        """Analyze secondary structure using DSSP."""
        system = self.get_system(system_name)
        protein = system.protein
        
        if protein is None:
            return
            
        try:
            frame_slice = self.get_frame_slice(system_name)
            
            # Run DSSP analysis
            dssp = secondary_structure.SecondaryStructure(protein.universe, protein.select_atoms("protein"))
            dssp.run(start=frame_slice[0], stop=frame_slice[-1]+1, step=frame_slice.step or 1)
            
            # Calculate secondary structure content
            ss_content = []
            times = self._get_time_array(system_name)
            
            for i, frame_idx in enumerate(frame_slice):
                frame_ss = dssp.results.timeseries[i]
                
                # Count secondary structure types
                helix_count = np.sum(np.isin(frame_ss, ['H', 'G', 'I']))  # α-helix, 3₁₀-helix, π-helix
                sheet_count = np.sum(np.isin(frame_ss, ['E', 'B']))       # β-sheet, β-bridge
                coil_count = np.sum(np.isin(frame_ss, ['T', 'S', 'C']))   # turn, bend, coil
                
                total_res = len(frame_ss)
                
                ss_content.append({
                    "helix_frac": helix_count / total_res if total_res > 0 else 0,
                    "sheet_frac": sheet_count / total_res if total_res > 0 else 0,
                    "coil_frac": coil_count / total_res if total_res > 0 else 0,
                    "helix_count": helix_count,
                    "sheet_count": sheet_count,
                    "coil_count": coil_count
                })
            
            # Create DataFrame
            df = pd.DataFrame(ss_content)
            df.index = times
            df.index.name = "time_ps"
            
            if "secondary_structure" not in self.results:
                self.results["secondary_structure"] = {}
            self.results["secondary_structure"][system_name] = df
            
            # Store raw DSSP results
            if "dssp_raw" not in self.results:
                self.results["dssp_raw"] = {}
            self.results["dssp_raw"][system_name] = dssp.results.timeseries
            
            logger.info(f"Secondary structure analysis completed for {system_name}")
            
        except Exception as e:
            logger.warning(f"Secondary structure analysis failed for {system_name}: {e}")
    
    def _analyze_radius_of_gyration(self, system_name: str) -> None:
        """Calculate radius of gyration."""
        system = self.get_system(system_name)
        protein = system.protein
        
        if protein is None:
            return
            
        frame_slice = self.get_frame_slice(system_name)
        
        rg_values = []
        rg_components = []
        
        for frame_idx in tqdm(frame_slice, desc=f"Radius of gyration {system_name}"):
            system.universe.trajectory[frame_idx]
            
            # Overall radius of gyration
            rg = protein.radius_of_gyration()
            rg_values.append(rg)
            
            # Components (Rg in each dimension)
            positions = protein.positions
            center_of_mass = protein.center_of_mass()
            
            # Calculate Rg components
            rel_pos = positions - center_of_mass
            masses = protein.masses
            total_mass = np.sum(masses)
            
            rg_x = np.sqrt(np.sum(masses * rel_pos[:, 0]**2) / total_mass)
            rg_y = np.sqrt(np.sum(masses * rel_pos[:, 1]**2) / total_mass)
            rg_z = np.sqrt(np.sum(masses * rel_pos[:, 2]**2) / total_mass)
            
            rg_components.append([rg_x, rg_y, rg_z])
        
        # Create DataFrame
        rg_components = np.array(rg_components)
        data = {
            "radius_of_gyration": np.array(rg_values),
            "rg_x": rg_components[:, 0],
            "rg_y": rg_components[:, 1], 
            "rg_z": rg_components[:, 2]
        }
        
        df = self._create_dataframe(data, system_name)
        
        if "radius_of_gyration" not in self.results:
            self.results["radius_of_gyration"] = {}
        self.results["radius_of_gyration"][system_name] = df
        
        logger.info(f"Radius of gyration analysis completed for {system_name}")
    
    def _analyze_end_to_end_distance(self, system_name: str) -> None:
        """Calculate end-to-end distance."""
        system = self.get_system(system_name)
        protein = system.protein
        
        if protein is None:
            return
            
        # Get N and C termini
        residues = protein.residues
        n_terminus = residues[0].atoms.select_atoms("name CA")[0]
        c_terminus = residues[-1].atoms.select_atoms("name CA")[0]
        
        frame_slice = self.get_frame_slice(system_name)
        
        distances = []
        
        for frame_idx in tqdm(frame_slice, desc=f"End-to-end distance {system_name}"):
            system.universe.trajectory[frame_idx]
            
            # Calculate distance
            distance = np.linalg.norm(c_terminus.position - n_terminus.position)
            distances.append(distance)
        
        # Create DataFrame
        data = {"end_to_end_distance": np.array(distances)}
        df = self._create_dataframe(data, system_name)
        
        if "end_to_end_distance" not in self.results:
            self.results["end_to_end_distance"] = {}
        self.results["end_to_end_distance"][system_name] = df
        
        logger.info(f"End-to-end distance analysis completed for {system_name}")
    
    def _analyze_contact_maps(self, system_name: str) -> None:
        """Calculate residue contact maps."""
        system = self.get_system(system_name)
        protein = system.protein
        
        if protein is None:
            return
            
        # Use CA atoms for contact analysis
        ca_atoms = protein.select_atoms("name CA")
        
        frame_slice = self.get_frame_slice(system_name)
        n_frames = len(frame_slice)
        n_residues = len(ca_atoms)
        
        # Initialize contact matrix
        contact_matrix = np.zeros((n_residues, n_residues))
        
        cutoff = 8.0  # Angstrom
        
        for frame_idx in tqdm(frame_slice, desc=f"Contact maps {system_name}"):
            system.universe.trajectory[frame_idx]
            
            # Calculate distance matrix
            positions = ca_atoms.positions
            dist_matrix = np.linalg.norm(
                positions[:, np.newaxis, :] - positions[np.newaxis, :, :], 
                axis=2
            )
            
            # Contacts (excluding neighbors)
            contact_frame = (dist_matrix <= cutoff) & (dist_matrix > 0)
            
            # Exclude immediate neighbors (|i-j| <= 3)
            for i in range(n_residues):
                for j in range(max(0, i-3), min(n_residues, i+4)):
                    contact_frame[i, j] = False
            
            contact_matrix += contact_frame
        
        # Normalize by number of frames
        contact_matrix = contact_matrix / n_frames
        
        if "contact_maps" not in self.results:
            self.results["contact_maps"] = {}
        self.results["contact_maps"][system_name] = contact_matrix
        
        logger.info(f"Contact map analysis completed for {system_name}")
    
    def get_structural_features(self, system_name: str) -> pd.DataFrame:
        """Get combined structural features for clustering.
        
        Args:
            system_name: Name of system
            
        Returns:
            DataFrame with structural features
        """
        features = {}
        
        # Radius of gyration
        if "radius_of_gyration" in self.results and system_name in self.results["radius_of_gyration"]:
            rg_data = self.results["radius_of_gyration"][system_name]
            features["radius_of_gyration"] = rg_data["radius_of_gyration"].values
        
        # End-to-end distance
        if "end_to_end_distance" in self.results and system_name in self.results["end_to_end_distance"]:
            ete_data = self.results["end_to_end_distance"][system_name]
            features["end_to_end_distance"] = ete_data["end_to_end_distance"].values
        
        # Secondary structure
        if "secondary_structure" in self.results and system_name in self.results["secondary_structure"]:
            ss_data = self.results["secondary_structure"][system_name]
            features["helix_frac"] = ss_data["helix_frac"].values
            features["sheet_frac"] = ss_data["sheet_frac"].values
            features["coil_frac"] = ss_data["coil_frac"].values
        
        if not features:
            raise ValueError(f"No structural features available for {system_name}")
        
        # Ensure all features have same length
        lengths = [len(v) for v in features.values()]
        if len(set(lengths)) > 1:
            min_length = min(lengths)
            features = {k: v[:min_length] for k, v in features.items()}
        
        # Get times for index
        times = None
        if "radius_of_gyration" in self.results and system_name in self.results["radius_of_gyration"]:
            times = self.results["radius_of_gyration"][system_name].index[:len(features[list(features.keys())[0]])]
        
        df = pd.DataFrame(features)
        if times is not None:
            df.index = times
            df.index.name = "time_ps"
        
        return df 
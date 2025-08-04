"""Protein-membrane orientation analysis based on K-Ras4B methodology."""

from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
from tqdm import tqdm
from loguru import logger
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

from ..core.base import BaseAnalyzer


class OrientationAnalyzer(BaseAnalyzer):
    """Analyzes protein orientation relative to membrane."""
    
    def __init__(self, config, effector_lobe_residues: Optional[str] = None, 
                 allosteric_lobe_residues: Optional[str] = None,
                 beta4_residues: Optional[str] = None, 
                 beta5_residues: Optional[str] = None):
        """Initialize orientation analyzer.
        
        Args:
            config: Analysis configuration
            effector_lobe_residues: Selection string for effector lobe residues
            allosteric_lobe_residues: Selection string for allosteric lobe residues  
            beta4_residues: Selection string for β4 strand residues
            beta5_residues: Selection string for β5 strand residues
        """
        super().__init__(config)
        
        # Default selections for K-Ras4B (can be customized)
        self.effector_lobe_residues = effector_lobe_residues or "resid 1-86"
        self.allosteric_lobe_residues = allosteric_lobe_residues or "resid 87-166"
        self.beta4_residues = beta4_residues or "resid 55-65"  # Typical β4 in Ras
        self.beta5_residues = beta5_residues or "resid 75-83"  # Typical β5 in Ras
        
    def run_analysis(self, system_names: Optional[List[str]] = None) -> None:
        """Run protein-membrane orientation analysis.
        
        Args:
            system_names: Names of systems to analyze. If None, analyze all.
        """
        if system_names is None:
            system_names = [name for name, sys in self._systems.items() 
                           if sys.has_protein and sys.has_membrane]
        else:
            system_names = [name for name in system_names 
                           if name in self._systems 
                           and self._systems[name].has_protein 
                           and self._systems[name].has_membrane]
        
        if not system_names:
            logger.warning("No systems with both protein and membrane found")
            return
            
        self._validate_systems(system_names)
        
        for system_name in system_names:
            logger.info(f"Analyzing protein-membrane orientation in system: {system_name}")
            
            # Calculate Dz and theta orientation parameters
            self._analyze_membrane_orientation(system_name)
            
        logger.info("Protein-membrane orientation analysis completed")
        
    def _analyze_membrane_orientation(self, system_name: str) -> None:
        """Analyze protein orientation relative to membrane.
        
        Args:
            system_name: Name of system to analyze
        """
        system = self.get_system(system_name)
        protein = system.protein
        membrane = system.membrane
        
        if protein is None or membrane is None:
            logger.warning(f"System {system_name} missing protein or membrane")
            return
            
        frame_slice = self.get_frame_slice(system_name)
        
        # Storage for results
        dz_values = []
        theta_values = []
        vx_vectors = []
        vy_vectors = []
        vz_vectors = []
        membrane_com_z = []
        effector_com_z = []
        
        for frame_idx in tqdm(frame_slice, desc=f"Orientation analysis {system_name}"):
            system.universe.trajectory[frame_idx]
            
            # Calculate Dz: distance between effector lobe COM and membrane COM in Z direction
            dz, eff_com_z, mem_com_z = self._calculate_dz(protein, membrane)
            dz_values.append(dz)
            effector_com_z.append(eff_com_z)
            membrane_com_z.append(mem_com_z)
            
            # Calculate theta: tilt angle relative to membrane
            theta, vx, vy, vz = self._calculate_theta(protein)
            theta_values.append(theta)
            vx_vectors.append(vx)
            vy_vectors.append(vy)
            vz_vectors.append(vz)
        
        # Create DataFrame
        data = {
            "dz": np.array(dz_values),
            "theta_degrees": np.array(theta_values),
            "effector_com_z": np.array(effector_com_z),
            "membrane_com_z": np.array(membrane_com_z),
            "vx_x": np.array([v[0] for v in vx_vectors]),
            "vx_y": np.array([v[1] for v in vx_vectors]),
            "vx_z": np.array([v[2] for v in vx_vectors]),
            "vy_x": np.array([v[0] for v in vy_vectors]),
            "vy_y": np.array([v[1] for v in vy_vectors]),
            "vy_z": np.array([v[2] for v in vy_vectors]),
            "vz_x": np.array([v[0] for v in vz_vectors]),
            "vz_y": np.array([v[1] for v in vz_vectors]),
            "vz_z": np.array([v[2] for v in vz_vectors]),
        }
        
        df = self._create_dataframe(data, system_name)
        
        if "orientation" not in self.results:
            self.results["orientation"] = {}
        self.results["orientation"][system_name] = df
        
        logger.info(f"Orientation analysis completed for {system_name}")
        
    def _calculate_dz(self, protein, membrane) -> Tuple[float, float, float]:
        """Calculate Dz distance between effector lobe and membrane COMs.
        
        Args:
            protein: Protein atom group
            membrane: Membrane atom group
            
        Returns:
            Tuple of (dz_distance, effector_com_z, membrane_com_z)
        """
        # Select effector lobe atoms
        try:
            effector_lobe = protein.select_atoms(self.effector_lobe_residues)
        except Exception as e:
            logger.warning(f"Could not select effector lobe: {e}, using entire protein")
            effector_lobe = protein
            
        # Calculate center of mass
        effector_com = effector_lobe.center_of_mass()
        membrane_com = membrane.center_of_mass()
        
        # Calculate Z-direction distance
        dz = effector_com[2] - membrane_com[2]
        
        return dz, effector_com[2], membrane_com[2]
    
    def _calculate_theta(self, protein) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray]:
        """Calculate theta tilt angle relative to membrane.
        
        Args:
            protein: Protein atom group
            
        Returns:
            Tuple of (theta_degrees, vx_vector, vy_vector, vz_vector)
        """
        # Define direction vectors as described in the methodology
        
        # Vx: vector connecting effector lobe COM to allosteric lobe COM  
        try:
            effector_lobe = protein.select_atoms(self.effector_lobe_residues)
            allosteric_lobe = protein.select_atoms(self.allosteric_lobe_residues)
            
            effector_com = effector_lobe.center_of_mass()
            allosteric_com = allosteric_lobe.center_of_mass()
            vx = allosteric_com - effector_com
        except Exception as e:
            logger.warning(f"Could not calculate Vx vector: {e}, using fallback")
            # Fallback: use first to last residue
            residues = protein.residues
            if len(residues) > 1:
                first_com = residues[0].atoms.center_of_mass()
                last_com = residues[-1].atoms.center_of_mass()
                vx = last_com - first_com
            else:
                vx = np.array([1.0, 0.0, 0.0])  # Default X direction
        
        # Vy: average orientation of β4 and β5 strands (main axis direction)
        try:
            beta4 = protein.select_atoms(self.beta4_residues)
            beta5 = protein.select_atoms(self.beta5_residues)
            
            # Calculate average direction of beta strands
            if len(beta4) > 0 and len(beta5) > 0:
                beta4_com = beta4.center_of_mass()
                beta5_com = beta5.center_of_mass()
                
                # Direction from beta4 to beta5
                beta_direction = beta5_com - beta4_com
                
                # For Vy, we want the average orientation - use normalized direction
                vy = beta_direction
            else:
                logger.warning("Could not find beta strands, using fallback for Vy")
                vy = np.array([0.0, 1.0, 0.0])  # Default Y direction
        except Exception as e:
            logger.warning(f"Could not calculate Vy vector: {e}, using fallback")
            vy = np.array([0.0, 1.0, 0.0])  # Default Y direction
            
        # Normalize vectors
        if np.linalg.norm(vx) > 0:
            vx = vx / np.linalg.norm(vx)
        if np.linalg.norm(vy) > 0:
            vy = vy / np.linalg.norm(vy)
            
        # Vz: cross product of Vx and Vy (principal axis)
        vz = np.cross(vx, vy)
        if np.linalg.norm(vz) > 0:
            vz = vz / np.linalg.norm(vz)
        else:
            vz = np.array([0.0, 0.0, 1.0])  # Default Z direction
        
        # Calculate theta: angle between membrane normal (Mz) and Vz
        membrane_normal = np.array([0.0, 0.0, 1.0])  # Z-axis is membrane normal
        
        # Calculate angle using dot product
        cos_theta = np.dot(vz, membrane_normal)
        cos_theta = np.clip(cos_theta, -1.0, 1.0)  # Ensure valid range
        theta_rad = np.arccos(cos_theta)
        theta_deg = np.degrees(theta_rad)
        
        return theta_deg, vx, vy, vz
    
    def get_orientation_summary(self, system_name: str) -> Dict[str, float]:
        """Get summary statistics for orientation analysis.
        
        Args:
            system_name: Name of system
            
        Returns:
            Dictionary with summary statistics
        """
        if "orientation" not in self.results or system_name not in self.results["orientation"]:
            raise ValueError(f"No orientation results found for system: {system_name}")
            
        df = self.results["orientation"][system_name]
        
        return {
            "mean_dz": float(df["dz"].mean()),
            "std_dz": float(df["dz"].std()),
            "mean_theta": float(df["theta_degrees"].mean()),
            "std_theta": float(df["theta_degrees"].std()),
            "min_dz": float(df["dz"].min()),
            "max_dz": float(df["dz"].max()),
            "min_theta": float(df["theta_degrees"].min()),
            "max_theta": float(df["theta_degrees"].max()),
        }
    
    def plot_orientation_timeseries(self, system_name: str, save_path: Optional[str] = None):
        """Plot orientation parameters over time.
        
        Args:
            system_name: Name of system
            save_path: Path to save plot
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            logger.error("matplotlib not available for plotting")
            return
            
        if "orientation" not in self.results or system_name not in self.results["orientation"]:
            logger.error(f"No orientation results found for system: {system_name}")
            return
            
        df = self.results["orientation"][system_name]
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
        
        # Plot Dz
        ax1.plot(df.index, df["dz"], 'b-', alpha=0.7)
        ax1.set_ylabel("Dz (Å)")
        ax1.set_title(f"Protein-Membrane Orientation Analysis: {system_name}")
        ax1.grid(True, alpha=0.3)
        
        # Plot theta
        ax2.plot(df.index, df["theta_degrees"], 'r-', alpha=0.7)
        ax2.set_xlabel("Time (ps)")
        ax2.set_ylabel("θ (degrees)")
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Plot saved to {save_path}")
        else:
            plt.show() 
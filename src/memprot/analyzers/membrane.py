"""Membrane lipid bilayer analysis."""

from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
from tqdm import tqdm
from loguru import logger
import MDAnalysis as mda
from MDAnalysis.analysis import leaflet
from MDAnalysis.lib.distances import distance_array

from ..core.base import BaseAnalyzer


class MembraneAnalyzer(BaseAnalyzer):
    """Analyzes membrane properties and dynamics."""
    
    def run_analysis(self, system_names: Optional[List[str]] = None) -> None:
        """Run membrane analysis.
        
        Args:
            system_names: Names of systems to analyze. If None, analyze all.
        """
        if system_names is None:
            system_names = [name for name, sys in self._systems.items() if sys.has_membrane]
        else:
            system_names = [name for name in system_names 
                           if name in self._systems and self._systems[name].has_membrane]
        
        if not system_names:
            logger.warning("No systems with membrane found")
            return
            
        self._validate_systems(system_names)
        
        for system_name in system_names:
            logger.info(f"Analyzing membrane in system: {system_name}")
            
            # Membrane thickness
            self._analyze_membrane_thickness(system_name)
            
            # Lipid composition
            self._analyze_lipid_composition(system_name)
            
            # Area per lipid
            self._analyze_area_per_lipid(system_name)
            
            # Membrane curvature
            self._analyze_membrane_curvature(system_name)
            
            # Lipid order parameters
            self._analyze_order_parameters(system_name)
            
            # Density profiles
            self._analyze_density_profiles(system_name)
            
        logger.info("Membrane analysis completed")
    
    def _analyze_membrane_thickness(self, system_name: str) -> None:
        """Analyze membrane thickness over time."""
        system = self.get_system(system_name)
        membrane = system.membrane
        
        if membrane is None:
            return
            
        frame_slice = self.get_frame_slice(system_name)
        
        thicknesses = []
        leaflet_separations = []
        
        for frame_idx in tqdm(frame_slice, desc=f"Membrane thickness {system_name}"):
            system.universe.trajectory[frame_idx]
            
            # Get membrane Z coordinates
            membrane_z = membrane.positions[:, 2]
            
            # Simple thickness calculation
            thickness = np.max(membrane_z) - np.min(membrane_z)
            thicknesses.append(thickness)
            
            # Try to identify leaflets for more accurate calculation
            try:
                # Get phosphorus atoms (common in phospholipids)
                phosphorus = membrane.select_atoms("name P*")
                if len(phosphorus) > 0:
                    p_z = phosphorus.positions[:, 2]
                    
                    # Separate into upper and lower leaflets
                    z_center = np.mean(p_z)
                    upper_leaflet = p_z[p_z > z_center]
                    lower_leaflet = p_z[p_z <= z_center]
                    
                    if len(upper_leaflet) > 0 and len(lower_leaflet) > 0:
                        leaflet_sep = np.mean(upper_leaflet) - np.mean(lower_leaflet)
                        leaflet_separations.append(leaflet_sep)
                    else:
                        leaflet_separations.append(thickness)
                else:
                    leaflet_separations.append(thickness)
            except:
                leaflet_separations.append(thickness)
        
        # Create DataFrame
        data = {
            "membrane_thickness": np.array(thicknesses),
            "leaflet_separation": np.array(leaflet_separations)
        }
        
        df = self._create_dataframe(data, system_name)
        
        if "membrane_thickness" not in self.results:
            self.results["membrane_thickness"] = {}
        self.results["membrane_thickness"][system_name] = df
        
        logger.info(f"Membrane thickness analysis completed for {system_name}")
    
    def _analyze_lipid_composition(self, system_name: str) -> None:
        """Analyze lipid composition."""
        system = self.get_system(system_name)
        membrane = system.membrane
        
        if membrane is None:
            return
            
        # Get unique residue names (lipid types)
        lipid_types = {}
        
        for residue in membrane.residues:
            resname = residue.resname
            if resname not in lipid_types:
                lipid_types[resname] = 0
            lipid_types[resname] += 1
        
        # Calculate mole fractions
        total_lipids = sum(lipid_types.values())
        mole_fractions = {lipid: count/total_lipids for lipid, count in lipid_types.items()}
        
        if "lipid_composition" not in self.results:
            self.results["lipid_composition"] = {}
        
        self.results["lipid_composition"][system_name] = {
            "counts": lipid_types,
            "mole_fractions": mole_fractions,
            "total_lipids": total_lipids
        }
        
        logger.info(f"Lipid composition analysis completed for {system_name}")
    
    def _analyze_area_per_lipid(self, system_name: str) -> None:
        """Analyze area per lipid."""
        system = self.get_system(system_name)
        membrane = system.membrane
        
        if membrane is None:
            return
            
        frame_slice = self.get_frame_slice(system_name)
        
        areas_per_lipid = []
        
        # Try to get phosphorus atoms for area calculation
        try:
            phosphorus = membrane.select_atoms("name P*")
            n_lipids = len(phosphorus)
            
            if n_lipids == 0:
                # Fallback: count residues
                n_lipids = len(membrane.residues)
                
        except:
            n_lipids = len(membrane.residues)
        
        for frame_idx in tqdm(frame_slice, desc=f"Area per lipid {system_name}"):
            system.universe.trajectory[frame_idx]
            
            # Get box dimensions
            box = system.universe.dimensions
            if box is not None and len(box) >= 3:
                # Assuming membrane in XY plane
                area = box[0] * box[1]  # A^2
                
                # Area per lipid (divide by 2 for bilayer)
                area_per_lipid = area / (n_lipids / 2) if n_lipids > 0 else 0
                areas_per_lipid.append(area_per_lipid)
            else:
                areas_per_lipid.append(np.nan)
        
        # Create DataFrame
        data = {"area_per_lipid": np.array(areas_per_lipid)}
        df = self._create_dataframe(data, system_name)
        
        if "area_per_lipid" not in self.results:
            self.results["area_per_lipid"] = {}
        self.results["area_per_lipid"][system_name] = df
        
        logger.info(f"Area per lipid analysis completed for {system_name}")
    
    def _analyze_membrane_curvature(self, system_name: str) -> None:
        """Analyze membrane curvature (simplified approach)."""
        system = self.get_system(system_name)
        membrane = system.membrane
        
        if membrane is None:
            return
            
        # Try to get phosphorus atoms for curvature analysis
        try:
            phosphorus = membrane.select_atoms("name P*")
            if len(phosphorus) == 0:
                # Fallback: use all membrane atoms
                phosphorus = membrane
        except:
            phosphorus = membrane
            
        frame_slice = self.get_frame_slice(system_name)
        
        curvatures = []
        undulations = []
        
        for frame_idx in tqdm(frame_slice, desc=f"Membrane curvature {system_name}"):
            system.universe.trajectory[frame_idx]
            
            positions = phosphorus.positions
            z_coords = positions[:, 2]
            
            # Calculate Z-coordinate variance as measure of undulation
            z_var = np.var(z_coords)
            undulations.append(z_var)
            
            # Simple curvature measure: standard deviation of Z
            z_std = np.std(z_coords)
            curvatures.append(z_std)
        
        # Create DataFrame
        data = {
            "membrane_curvature": np.array(curvatures),
            "membrane_undulation": np.array(undulations)
        }
        
        df = self._create_dataframe(data, system_name)
        
        if "membrane_curvature" not in self.results:
            self.results["membrane_curvature"] = {}
        self.results["membrane_curvature"][system_name] = df
        
        logger.info(f"Membrane curvature analysis completed for {system_name}")
    
    def _analyze_order_parameters(self, system_name: str) -> None:
        """Analyze lipid order parameters (simplified)."""
        system = self.get_system(system_name)
        membrane = system.membrane
        
        if membrane is None:
            return
            
        frame_slice = self.get_frame_slice(system_name)
        
        # Try to find carbon chains for order parameter calculation
        try:
            # Look for carbon atoms in acyl chains
            carbon_chains = membrane.select_atoms("name C* and not name C1*")
            
            if len(carbon_chains) == 0:
                logger.warning(f"No carbon chain atoms found in {system_name}")
                return
                
        except:
            logger.warning(f"Could not select carbon atoms in {system_name}")
            return
        
        order_params = []
        
        for frame_idx in tqdm(frame_slice, desc=f"Order parameters {system_name}"):
            system.universe.trajectory[frame_idx]
            
            # Simple order parameter: alignment with membrane normal (Z-axis)
            positions = carbon_chains.positions
            
            # Calculate vectors between consecutive carbons
            if len(positions) > 1:
                vectors = positions[1:] - positions[:-1]
                
                # Calculate angle with Z-axis
                z_unit = np.array([0, 0, 1])
                cos_angles = []
                
                for vec in vectors:
                    if np.linalg.norm(vec) > 0:
                        cos_angle = np.dot(vec, z_unit) / np.linalg.norm(vec)
                        cos_angles.append(cos_angle)
                
                if cos_angles:
                    # Order parameter S = <3cos²θ - 1>/2
                    cos_angles = np.array(cos_angles)
                    order_param = np.mean(1.5 * cos_angles**2 - 0.5)
                    order_params.append(order_param)
                else:
                    order_params.append(0)
            else:
                order_params.append(0)
        
        # Create DataFrame
        data = {"order_parameter": np.array(order_params)}
        df = self._create_dataframe(data, system_name)
        
        if "order_parameters" not in self.results:
            self.results["order_parameters"] = {}
        self.results["order_parameters"][system_name] = df
        
        logger.info(f"Order parameter analysis completed for {system_name}")
    
    def _analyze_density_profiles(self, system_name: str) -> None:
        """Analyze membrane density profiles along Z-axis."""
        system = self.get_system(system_name)
        membrane = system.membrane
        
        if membrane is None:
            return
            
        frame_slice = self.get_frame_slice(system_name)
        
        # Define Z-bins
        n_bins = 50
        
        # Collect Z coordinates from all frames
        all_z_coords = []
        z_min_global = float('inf')
        z_max_global = float('-inf')
        
        # First pass: find global Z range
        for frame_idx in frame_slice:
            system.universe.trajectory[frame_idx]
            z_coords = membrane.positions[:, 2]
            z_min_global = min(z_min_global, np.min(z_coords))
            z_max_global = max(z_max_global, np.max(z_coords))
        
        # Create bins
        z_bins = np.linspace(z_min_global, z_max_global, n_bins + 1)
        z_centers = (z_bins[:-1] + z_bins[1:]) / 2
        
        density_profiles = []
        
        # Second pass: calculate density profiles
        for frame_idx in tqdm(frame_slice, desc=f"Density profiles {system_name}"):
            system.universe.trajectory[frame_idx]
            
            z_coords = membrane.positions[:, 2]
            
            # Calculate histogram
            counts, _ = np.histogram(z_coords, bins=z_bins)
            
            # Normalize by bin volume (bin_width * box_area)
            box = system.universe.dimensions
            if box is not None and len(box) >= 3:
                bin_width = z_bins[1] - z_bins[0]
                box_area = box[0] * box[1]
                bin_volume = bin_width * box_area
                
                density = counts / bin_volume  # atoms/Å³
            else:
                density = counts
            
            density_profiles.append(density)
        
        # Average density profile
        density_profiles = np.array(density_profiles)
        avg_density = np.mean(density_profiles, axis=0)
        std_density = np.std(density_profiles, axis=0)
        
        # Store results
        if "density_profiles" not in self.results:
            self.results["density_profiles"] = {}
        
        self.results["density_profiles"][system_name] = {
            "z_centers": z_centers,
            "avg_density": avg_density,
            "std_density": std_density,
            "all_profiles": density_profiles
        }
        
        logger.info(f"Density profile analysis completed for {system_name}")
    
    def compare_membrane_properties(self, system1: str, system2: str) -> Dict[str, pd.DataFrame]:
        """Compare membrane properties between two systems.
        
        Args:
            system1: First system name
            system2: Second system name
            
        Returns:
            Dictionary with comparison DataFrames
        """
        comparisons = {}
        
        # Compare thickness
        if ("membrane_thickness" in self.results and 
            system1 in self.results["membrane_thickness"] and 
            system2 in self.results["membrane_thickness"]):
            
            df1 = self.results["membrane_thickness"][system1]
            df2 = self.results["membrane_thickness"][system2]
            
            comparison = pd.DataFrame({
                f"{system1}_thickness": df1["membrane_thickness"],
                f"{system2}_thickness": df2["membrane_thickness"],
                f"{system1}_leaflet_sep": df1["leaflet_separation"],
                f"{system2}_leaflet_sep": df2["leaflet_separation"],
            })
            
            comparisons["thickness"] = comparison
        
        # Compare area per lipid
        if ("area_per_lipid" in self.results and 
            system1 in self.results["area_per_lipid"] and 
            system2 in self.results["area_per_lipid"]):
            
            df1 = self.results["area_per_lipid"][system1]
            df2 = self.results["area_per_lipid"][system2]
            
            comparison = pd.DataFrame({
                f"{system1}_area": df1["area_per_lipid"],
                f"{system2}_area": df2["area_per_lipid"],
            })
            
            comparisons["area_per_lipid"] = comparison
        
        # Compare curvature
        if ("membrane_curvature" in self.results and 
            system1 in self.results["membrane_curvature"] and 
            system2 in self.results["membrane_curvature"]):
            
            df1 = self.results["membrane_curvature"][system1]
            df2 = self.results["membrane_curvature"][system2]
            
            comparison = pd.DataFrame({
                f"{system1}_curvature": df1["membrane_curvature"],
                f"{system2}_curvature": df2["membrane_curvature"],
                f"{system1}_undulation": df1["membrane_undulation"],
                f"{system2}_undulation": df2["membrane_undulation"],
            })
            
            comparisons["curvature"] = comparison
        
        return comparisons
    
    def get_membrane_influence_metrics(self, membrane_only_system: str, protein_systems: List[str]) -> Dict[str, Any]:
        """Calculate membrane influence metrics.
        
        Args:
            membrane_only_system: System with membrane only
            protein_systems: Systems with protein
            
        Returns:
            Dictionary with influence metrics
        """
        metrics = {}
        
        if membrane_only_system not in self._systems:
            logger.error(f"Membrane-only system {membrane_only_system} not found")
            return metrics
        
        # Get membrane-only baseline
        baseline_results = {}
        for analysis_type in ["membrane_thickness", "area_per_lipid", "membrane_curvature", "order_parameters"]:
            if analysis_type in self.results and membrane_only_system in self.results[analysis_type]:
                baseline_results[analysis_type] = self.results[analysis_type][membrane_only_system]
        
        # Compare with protein systems
        for protein_system in protein_systems:
            if protein_system not in self._systems:
                continue
                
            system_metrics = {}
            
            for analysis_type, baseline_df in baseline_results.items():
                if analysis_type in self.results and protein_system in self.results[analysis_type]:
                    protein_df = self.results[analysis_type][protein_system]
                    
                    # Calculate differences
                    if analysis_type == "membrane_thickness":
                        baseline_mean = baseline_df["membrane_thickness"].mean()
                        protein_mean = protein_df["membrane_thickness"].mean()
                        thickness_change = protein_mean - baseline_mean
                        thickness_change_pct = (thickness_change / baseline_mean) * 100
                        
                        system_metrics["thickness_change"] = thickness_change
                        system_metrics["thickness_change_pct"] = thickness_change_pct
                        
                    elif analysis_type == "area_per_lipid":
                        baseline_mean = baseline_df["area_per_lipid"].mean()
                        protein_mean = protein_df["area_per_lipid"].mean()
                        area_change = protein_mean - baseline_mean
                        area_change_pct = (area_change / baseline_mean) * 100
                        
                        system_metrics["area_change"] = area_change
                        system_metrics["area_change_pct"] = area_change_pct
                        
                    elif analysis_type == "membrane_curvature":
                        baseline_mean = baseline_df["membrane_curvature"].mean()
                        protein_mean = protein_df["membrane_curvature"].mean()
                        curvature_change = protein_mean - baseline_mean
                        curvature_change_pct = (curvature_change / baseline_mean) * 100 if baseline_mean > 0 else 0
                        
                        system_metrics["curvature_change"] = curvature_change
                        system_metrics["curvature_change_pct"] = curvature_change_pct
            
            metrics[protein_system] = system_metrics
        
        return metrics 
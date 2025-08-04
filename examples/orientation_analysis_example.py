#!/usr/bin/env python3
"""
Example script demonstrating K-Ras4B membrane orientation analysis.

This script shows how to use the OrientationAnalyzer to calculate:
1. Dz: Distance between effector lobe COM and membrane COM in Z direction
2. θ (theta): Tilt angle of protein core domain relative to membrane

Based on the methodology described in the research article.
"""

import sys
from pathlib import Path
from memprot.core.config import AnalysisConfig
from memprot.core.system import SimulationSystem
from memprot.analyzers import OrientationAnalyzer

def main():
    """Run K-Ras4B orientation analysis example."""
    
    # Configure analysis parameters
    config = AnalysisConfig(
        start_time=100000,  # Start after 100 ns (100,000 ps) as mentioned in article
        end_time=None,      # Use all available frames after start_time
        dt=10,              # Analysis frequency (ps)
        output_dir="orientation_results"
    )
    
    # Initialize the orientation analyzer with custom residue selections
    # These selections should be adjusted based on your specific protein structure
    analyzer = OrientationAnalyzer(
        config=config,
        effector_lobe_residues="resid 1-86",    # K-Ras4B effector lobe
        allosteric_lobe_residues="resid 87-166", # K-Ras4B allosteric lobe  
        beta4_residues="resid 55-65",            # β4 strand residues
        beta5_residues="resid 75-83"             # β5 strand residues
    )
    
    # Example: Load your simulation system
    # Replace these paths with your actual topology and trajectory files
    try:
        system = SimulationSystem(
            name="kras4b_membrane",
            topology_file="path/to/your/topology.pdb",  # or .gro, .psf, etc.
            trajectory_file="path/to/your/trajectory.xtc",  # or .dcd, .trr, etc.
            protein_selection="protein",
            membrane_selection="resname POPC POPE POPS or name P*",  # Adjust as needed
        )
        
        # Add system to analyzer
        analyzer.add_system(system)
        
        # Run the orientation analysis
        print("Running K-Ras4B membrane orientation analysis...")
        analyzer.run_analysis()
        
        # Get summary statistics
        summary = analyzer.get_orientation_summary("kras4b_membrane")
        print("\nOrientation Analysis Summary:")
        print(f"Mean Dz: {summary['mean_dz']:.2f} ± {summary['std_dz']:.2f} Å")
        print(f"Mean θ: {summary['mean_theta']:.2f} ± {summary['std_theta']:.2f}°")
        print(f"Dz range: {summary['min_dz']:.2f} to {summary['max_dz']:.2f} Å")
        print(f"θ range: {summary['min_theta']:.2f} to {summary['max_theta']:.2f}°")
        
        # Plot time series
        analyzer.plot_orientation_timeseries("kras4b_membrane", "orientation_timeseries.png")
        
        # Save results
        analyzer.save_results()
        
        # Access detailed results
        results_df = analyzer.results["orientation"]["kras4b_membrane"]
        print(f"\nDetailed results saved with {len(results_df)} time points")
        print("Available columns:", list(results_df.columns))
        
        # Example: Calculate correlation between Dz and theta
        correlation = results_df["dz"].corr(results_df["theta_degrees"])
        print(f"Correlation between Dz and θ: {correlation:.3f}")
        
    except Exception as e:
        print(f"Error loading simulation system: {e}")
        print("\nTo use this example:")
        print("1. Replace topology_file and trajectory_file with your actual file paths")
        print("2. Adjust protein_selection and membrane_selection for your system")
        print("3. Modify residue selections to match your protein structure")
        

def analyze_custom_protein():
    """Example for analyzing other proteins besides K-Ras4B."""
    
    config = AnalysisConfig(
        start_time=50000,  # Start after 50 ns
        dt=5,              # More frequent sampling
        output_dir="custom_orientation_results"
    )
    
    # For other proteins, you'll need to identify:
    # 1. Effector/allosteric domains or N/C-terminal regions
    # 2. Key secondary structure elements for orientation vector
    analyzer = OrientationAnalyzer(
        config=config,
        effector_lobe_residues="resid 1-100",     # Adjust for your protein
        allosteric_lobe_residues="resid 101-200", # Adjust for your protein
        beta4_residues="resid 80-90",             # Key beta strand 1
        beta5_residues="resid 120-130"            # Key beta strand 2
    )
    
    print("Custom protein orientation analysis configured.")
    print("Adjust residue selections based on your protein structure.")


if __name__ == "__main__":
    main()
    print("\n" + "="*60)
    analyze_custom_protein() 
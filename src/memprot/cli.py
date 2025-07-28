"""Command-line interface for MD analysis pipeline."""

import click
from pathlib import Path
from typing import List, Optional
from loguru import logger

from .core import AnalysisConfig, SimulationSystem, create_example_config
from .analyzers import ProteinAnalyzer, MembraneAnalyzer, InteractionAnalyzer, ClusterAnalyzer
from .visualization import PlotManager
from .utils import setup_logging


@click.group()
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose logging')
@click.option('--log-file', type=click.Path(), help='Log file path')
def main(verbose: bool, log_file: Optional[str]) -> None:
    """Membrane protein analysis pipeline."""
    level = "DEBUG" if verbose else "INFO"
    setup_logging(level=level, log_file=log_file)
    logger.info("Starting membrane protein analysis pipeline")


@main.command()
@click.option('--config', '-c', type=click.Path(exists=True), required=True,
              help='Configuration file path')
@click.option('--systems', multiple=True, help='Specific systems to analyze')
@click.option('--skip-protein', is_flag=True, help='Skip protein analysis')
@click.option('--skip-membrane', is_flag=True, help='Skip membrane analysis')
@click.option('--skip-interactions', is_flag=True, help='Skip interaction analysis')
@click.option('--skip-clustering', is_flag=True, help='Skip clustering analysis')
@click.option('--skip-plots', is_flag=True, help='Skip generating plots')
def analyze(
    config: str,
    systems: tuple,
    skip_protein: bool,
    skip_membrane: bool,
    skip_interactions: bool,
    skip_clustering: bool,
    skip_plots: bool
) -> None:
    """Run comprehensive MD analysis."""
    try:
        # Load configuration
        config_obj = AnalysisConfig.from_yaml(config)
        config_obj.validate()
        
        # Filter systems if specified
        system_names = list(systems) if systems else list(config_obj.systems.keys())
        
        # Load simulation systems
        simulation_systems = {}
        for name in system_names:
            if name not in config_obj.systems:
                logger.error(f"System '{name}' not found in configuration")
                continue
                
            system = SimulationSystem(config_obj.systems[name])
            system.load()
            simulation_systems[name] = system
            logger.info(f"Loaded system: {name}")
        
        if not simulation_systems:
            logger.error("No valid systems loaded")
            return
        
        # Initialize analyzers
        analyzers = {}
        all_results = {}
        
        # Protein analysis
        if not skip_protein:
            logger.info("Running protein analysis...")
            protein_analyzer = ProteinAnalyzer(config_obj)
            for system in simulation_systems.values():
                protein_analyzer.add_system(system)
            
            protein_analyzer.run_analysis(system_names)
            protein_analyzer.save_results()
            analyzers['protein'] = protein_analyzer
            all_results.update(protein_analyzer.results)
        
        # Membrane analysis
        if not skip_membrane:
            logger.info("Running membrane analysis...")
            membrane_analyzer = MembraneAnalyzer(config_obj)
            for system in simulation_systems.values():
                membrane_analyzer.add_system(system)
            
            membrane_analyzer.run_analysis(system_names)
            membrane_analyzer.save_results()
            analyzers['membrane'] = membrane_analyzer
            all_results.update(membrane_analyzer.results)
        
        # Interaction analysis
        if not skip_interactions:
            logger.info("Running interaction analysis...")
            interaction_analyzer = InteractionAnalyzer(config_obj)
            for system in simulation_systems.values():
                interaction_analyzer.add_system(system)
            
            interaction_analyzer.run_analysis(system_names)
            interaction_analyzer.save_results()
            analyzers['interactions'] = interaction_analyzer
            all_results.update(interaction_analyzer.results)
        
        # Clustering analysis
        if not skip_clustering and 'protein' in analyzers:
            logger.info("Running clustering analysis...")
            cluster_analyzer = ClusterAnalyzer(config_obj)
            for system in simulation_systems.values():
                cluster_analyzer.add_system(system)
            
            # Pass protein results to clustering analyzer
            cluster_analyzer._protein_analyzer_results = analyzers['protein'].results
            cluster_analyzer.run_analysis(system_names)
            cluster_analyzer.save_results()
            analyzers['clustering'] = cluster_analyzer
            all_results.update(cluster_analyzer.results)
        
        # Generate plots
        if not skip_plots:
            logger.info("Generating plots...")
            plot_manager = PlotManager(config_obj)
            
            # Plot protein results
            if 'protein' in analyzers:
                if 'rmsd' in analyzers['protein'].results:
                    plot_manager.plot_protein_rmsd(analyzers['protein'].results['rmsd'])
                
                if 'radius_of_gyration' in analyzers['protein'].results:
                    plot_manager.plot_radius_of_gyration(analyzers['protein'].results['radius_of_gyration'])
                
                if 'secondary_structure' in analyzers['protein'].results:
                    plot_manager.plot_secondary_structure(analyzers['protein'].results['secondary_structure'])
                
                if 'contact_maps' in analyzers['protein'].results:
                    for system_name, contact_matrix in analyzers['protein'].results['contact_maps'].items():
                        plot_manager.plot_contact_map(contact_matrix, system_name)
            
            # Plot interaction results
            if 'interactions' in analyzers:
                if 'pip2_interactions' in analyzers['interactions'].results:
                    plot_manager.plot_pip2_interactions(analyzers['interactions'].results['pip2_interactions'])
            
            # Plot membrane results
            if 'membrane' in analyzers:
                thickness_results = analyzers['membrane'].results.get('membrane_thickness', {})
                area_results = analyzers['membrane'].results.get('area_per_lipid', {})
                
                if thickness_results:
                    plot_manager.plot_membrane_properties(thickness_results, area_results)
                
                if 'density_profiles' in analyzers['membrane'].results:
                    plot_manager.plot_density_profile(analyzers['membrane'].results['density_profiles'])
            
            # Plot clustering results
            if 'clustering' in analyzers and 'clustering_results' in analyzers['clustering'].results:
                plot_manager.plot_clustering_results(analyzers['clustering'].results['clustering_results'])
            
            # Create interactive dashboard
            dashboard = plot_manager.create_interactive_dashboard(all_results)
            plot_manager.save_interactive_dashboard(dashboard)
        
        # Generate summary report
        _generate_summary_report(config_obj, analyzers, system_names)
        
        logger.info("Analysis completed successfully!")
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        raise


@main.command()
@click.option('--output', '-o', type=click.Path(), default='config.yaml',
              help='Output configuration file path')
def create_config(output: str) -> None:
    """Create example configuration file."""
    config = create_example_config()
    config.to_yaml(output)
    logger.info(f"Example configuration created: {output}")


@main.command()
@click.option('--config', '-c', type=click.Path(exists=True), required=True,
              help='Configuration file path')
def validate_config(config: str) -> None:
    """Validate configuration file."""
    try:
        config_obj = AnalysisConfig.from_yaml(config)
        config_obj.validate()
        logger.info("Configuration is valid!")
    except Exception as e:
        logger.error(f"Configuration validation failed: {e}")


@main.command()
@click.option('--config', '-c', type=click.Path(exists=True), required=True,
              help='Configuration file path')
@click.option('--system1', required=True, help='First system name')
@click.option('--system2', required=True, help='Second system name')
def compare(config: str, system1: str, system2: str) -> None:
    """Compare analysis results between two systems."""
    try:
        config_obj = AnalysisConfig.from_yaml(config)
        
        # Load and analyze both systems
        systems = [system1, system2]
        simulation_systems = {}
        
        for name in systems:
            if name not in config_obj.systems:
                logger.error(f"System '{name}' not found in configuration")
                return
                
            system = SimulationSystem(config_obj.systems[name])
            system.load()
            simulation_systems[name] = system
        
        # Run quick analysis for comparison
        analyzers = {}
        
        # Protein analysis
        protein_analyzer = ProteinAnalyzer(config_obj)
        for system in simulation_systems.values():
            protein_analyzer.add_system(system)
        protein_analyzer.run_analysis(systems)
        analyzers['protein'] = protein_analyzer
        
        # Interaction analysis
        interaction_analyzer = InteractionAnalyzer(config_obj)
        for system in simulation_systems.values():
            interaction_analyzer.add_system(system)
        interaction_analyzer.run_analysis(systems)
        analyzers['interactions'] = interaction_analyzer
        
        # Generate comparison plots
        plot_manager = PlotManager(config_obj)
        
        # Compare PIP2 interactions
        if 'pip2_interactions' in interaction_analyzer.results:
            comparisons = interaction_analyzer.compare_pip2_interactions(system1, system2)
            for comparison_type, df in comparisons.items():
                logger.info(f"\n{comparison_type.upper()} COMPARISON:")
                logger.info(f"Mean {system1}: {df.iloc[:, 0].mean():.3f}")
                logger.info(f"Mean {system2}: {df.iloc[:, 1].mean():.3f}")
        
        # Plot comparison
        comparison_data = {
            system1: analyzers['protein'].get_summary_statistics(),
            system2: analyzers['protein'].get_summary_statistics()
        }
        plot_manager.plot_comparison_summary(comparison_data)
        
        logger.info(f"Comparison between {system1} and {system2} completed!")
        
    except Exception as e:
        logger.error(f"Comparison failed: {e}")
        raise


def _generate_summary_report(
    config: AnalysisConfig, 
    analyzers: dict, 
    system_names: List[str]
) -> None:
    """Generate summary report."""
    report_path = config.output_dir / "analysis_summary.txt"
    
    with open(report_path, 'w') as f:
        f.write("MD ANALYSIS SUMMARY REPORT\n")
        f.write("=" * 50 + "\n\n")
        
        f.write(f"Configuration: {config}\n")
        f.write(f"Systems analyzed: {', '.join(system_names)}\n\n")
        
        # Summary statistics for each analyzer
        for analyzer_name, analyzer in analyzers.items():
            f.write(f"{analyzer_name.upper()} ANALYSIS\n")
            f.write("-" * 30 + "\n")
            
            try:
                summary = analyzer.get_summary_statistics()
                for key, stats in summary.items():
                    f.write(f"\n{key}:\n")
                    if isinstance(stats, dict):
                        for stat_name, value in stats.items():
                            if isinstance(value, dict):
                                for subkey, subvalue in value.items():
                                    f.write(f"  {stat_name} {subkey}: {subvalue:.3f}\n")
                            else:
                                f.write(f"  {stat_name}: {value:.3f}\n")
            except Exception as e:
                f.write(f"Error generating summary: {e}\n")
            
            f.write("\n")
    
    logger.info(f"Summary report saved to {report_path}")


if __name__ == '__main__':
    main() 
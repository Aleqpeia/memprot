"""Plotting and visualization tools."""

from typing import Dict, List, Optional, Tuple, Any, Union
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from pathlib import Path
from loguru import logger

from ..core.config import AnalysisConfig


class PlotManager:
    """Manages plotting and visualization of analysis results."""
    
    def __init__(self, config: AnalysisConfig):
        """Initialize plot manager.
        
        Args:
            config: Analysis configuration
        """
        self.config = config
        self.output_dir = config.output_dir / "plots"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set style
        sns.set_style("whitegrid")
        plt.rcParams.update({'font.size': 12})
    
    def plot_protein_rmsd(self, results: Dict[str, pd.DataFrame], save: bool = True) -> plt.Figure:
        """Plot protein RMSD trajectories.
        
        Args:
            results: Dictionary of system_name -> RMSD DataFrame
            save: Whether to save the plot
            
        Returns:
            Matplotlib figure
        """
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
        
        colors = plt.cm.tab10(np.linspace(0, 1, len(results)))
        
        for i, (system_name, df) in enumerate(results.items()):
            time = df.index / 1000  # Convert to ns
            
            # CA RMSD
            ax1.plot(time, df['rmsd_ca'], label=system_name, color=colors[i], alpha=0.8)
            
            # Backbone RMSD
            ax2.plot(time, df['rmsd_backbone'], label=system_name, color=colors[i], alpha=0.8)
        
        ax1.set_ylabel('CA RMSD (Å)')
        ax1.legend()
        ax1.set_title('C-alpha RMSD')
        
        ax2.set_xlabel('Time (ns)')
        ax2.set_ylabel('Backbone RMSD (Å)')
        ax2.legend()
        ax2.set_title('Backbone RMSD')
        
        plt.tight_layout()
        
        if save and self.config.save_plots:
            save_path = self.output_dir / f"protein_rmsd.{self.config.plot_format}"
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"RMSD plot saved to {save_path}")
        
        return fig
    
    def plot_radius_of_gyration(self, results: Dict[str, pd.DataFrame], save: bool = True) -> plt.Figure:
        """Plot radius of gyration trajectories.
        
        Args:
            results: Dictionary of system_name -> Rg DataFrame
            save: Whether to save the plot
            
        Returns:
            Matplotlib figure
        """
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        axes = axes.flatten()
        
        colors = plt.cm.tab10(np.linspace(0, 1, len(results)))
        
        for i, (system_name, df) in enumerate(results.items()):
            time = df.index / 1000  # Convert to ns
            
            # Overall Rg
            axes[0].plot(time, df['radius_of_gyration'], label=system_name, color=colors[i], alpha=0.8)
            
            # Rg components
            axes[1].plot(time, df['rg_x'], label=f"{system_name} X", color=colors[i], alpha=0.8)
            axes[2].plot(time, df['rg_y'], label=f"{system_name} Y", color=colors[i], alpha=0.8)
            axes[3].plot(time, df['rg_z'], label=f"{system_name} Z", color=colors[i], alpha=0.8)
        
        axes[0].set_ylabel('Radius of Gyration (Å)')
        axes[0].set_title('Overall Radius of Gyration')
        axes[0].legend()
        
        for i, comp in enumerate(['X', 'Y', 'Z'], 1):
            axes[i].set_ylabel(f'Rg {comp} (Å)')
            axes[i].set_title(f'Radius of Gyration - {comp} Component')
            axes[i].legend()
            if i >= 2:
                axes[i].set_xlabel('Time (ns)')
        
        plt.tight_layout()
        
        if save and self.config.save_plots:
            save_path = self.output_dir / f"radius_of_gyration.{self.config.plot_format}"
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Radius of gyration plot saved to {save_path}")
        
        return fig
    
    def plot_secondary_structure(self, results: Dict[str, pd.DataFrame], save: bool = True) -> plt.Figure:
        """Plot secondary structure evolution.
        
        Args:
            results: Dictionary of system_name -> SS DataFrame
            save: Whether to save the plot
            
        Returns:
            Matplotlib figure
        """
        n_systems = len(results)
        fig, axes = plt.subplots(n_systems, 1, figsize=(12, 4*n_systems))
        
        if n_systems == 1:
            axes = [axes]
        
        for i, (system_name, df) in enumerate(results.items()):
            time = df.index / 1000  # Convert to ns
            
            # Stacked area plot
            axes[i].fill_between(time, 0, df['helix_frac'], alpha=0.7, label='Helix', color='red')
            axes[i].fill_between(time, df['helix_frac'], 
                               df['helix_frac'] + df['sheet_frac'], 
                               alpha=0.7, label='Sheet', color='blue')
            axes[i].fill_between(time, df['helix_frac'] + df['sheet_frac'], 
                               1.0, alpha=0.7, label='Coil', color='green')
            
            axes[i].set_ylabel('Fraction')
            axes[i].set_title(f'Secondary Structure - {system_name}')
            axes[i].legend()
            axes[i].set_ylim(0, 1)
            
            if i == n_systems - 1:
                axes[i].set_xlabel('Time (ns)')
        
        plt.tight_layout()
        
        if save and self.config.save_plots:
            save_path = self.output_dir / f"secondary_structure.{self.config.plot_format}"
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Secondary structure plot saved to {save_path}")
        
        return fig
    
    def plot_pip2_interactions(self, results: Dict[str, pd.DataFrame], save: bool = True) -> plt.Figure:
        """Plot PIP2 interaction analysis.
        
        Args:
            results: Dictionary of system_name -> PIP2 interaction DataFrame
            save: Whether to save the plot
            
        Returns:
            Matplotlib figure
        """
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        
        colors = plt.cm.tab10(np.linspace(0, 1, len(results)))
        
        for i, (system_name, df) in enumerate(results.items()):
            time = df.index / 1000  # Convert to ns
            
            # PIP2 contacts
            axes[0, 0].plot(time, df['pip2_contacts'], label=system_name, color=colors[i], alpha=0.8)
            
            # Average distance
            axes[0, 1].plot(time, df['avg_pip2_distance'], label=system_name, color=colors[i], alpha=0.8)
            
            # Closest distance
            axes[1, 0].plot(time, df['closest_pip2_distance'], label=system_name, color=colors[i], alpha=0.8)
            
            # Distribution of contacts
            axes[1, 1].hist(df['pip2_contacts'], bins=20, alpha=0.6, label=system_name, color=colors[i])
        
        axes[0, 0].set_ylabel('Number of Contacts')
        axes[0, 0].set_title('PIP2 Contacts')
        axes[0, 0].legend()
        
        axes[0, 1].set_ylabel('Distance (Å)')
        axes[0, 1].set_title('Average PIP2 Distance')
        axes[0, 1].legend()
        
        axes[1, 0].set_xlabel('Time (ns)')
        axes[1, 0].set_ylabel('Distance (Å)')
        axes[1, 0].set_title('Closest PIP2 Distance')
        axes[1, 0].legend()
        
        axes[1, 1].set_xlabel('Number of Contacts')
        axes[1, 1].set_ylabel('Frequency')
        axes[1, 1].set_title('Contact Distribution')
        axes[1, 1].legend()
        
        plt.tight_layout()
        
        if save and self.config.save_plots:
            save_path = self.output_dir / f"pip2_interactions.{self.config.plot_format}"
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"PIP2 interactions plot saved to {save_path}")
        
        return fig
    
    def plot_membrane_properties(self, thickness_results: Dict[str, pd.DataFrame], 
                                area_results: Optional[Dict[str, pd.DataFrame]] = None,
                                save: bool = True) -> plt.Figure:
        """Plot membrane properties.
        
        Args:
            thickness_results: Dictionary of system_name -> thickness DataFrame
            area_results: Dictionary of system_name -> area per lipid DataFrame
            save: Whether to save the plot
            
        Returns:
            Matplotlib figure
        """
        n_plots = 2 if area_results else 1
        fig, axes = plt.subplots(n_plots, 1, figsize=(12, 6*n_plots))
        
        if n_plots == 1:
            axes = [axes]
        
        colors = plt.cm.tab10(np.linspace(0, 1, len(thickness_results)))
        
        # Membrane thickness
        for i, (system_name, df) in enumerate(thickness_results.items()):
            time = df.index / 1000  # Convert to ns
            axes[0].plot(time, df['membrane_thickness'], label=system_name, color=colors[i], alpha=0.8)
        
        axes[0].set_ylabel('Thickness (Å)')
        axes[0].set_title('Membrane Thickness')
        axes[0].legend()
        
        # Area per lipid
        if area_results:
            for i, (system_name, df) in enumerate(area_results.items()):
                time = df.index / 1000  # Convert to ns
                axes[1].plot(time, df['area_per_lipid'], label=system_name, color=colors[i], alpha=0.8)
            
            axes[1].set_xlabel('Time (ns)')
            axes[1].set_ylabel('Area per Lipid (Å²)')
            axes[1].set_title('Area per Lipid')
            axes[1].legend()
        else:
            axes[0].set_xlabel('Time (ns)')
        
        plt.tight_layout()
        
        if save and self.config.save_plots:
            save_path = self.output_dir / f"membrane_properties.{self.config.plot_format}"
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Membrane properties plot saved to {save_path}")
        
        return fig
    
    def plot_clustering_results(self, results: Dict[str, Dict], method: str = "kmeans", save: bool = True) -> plt.Figure:
        """Plot clustering analysis results.
        
        Args:
            results: Clustering results from ClusterAnalyzer
            method: Clustering method to plot
            save: Whether to save the plot
            
        Returns:
            Matplotlib figure
        """
        n_systems = len(results)
        fig, axes = plt.subplots(2, n_systems, figsize=(5*n_systems, 10))
        
        if n_systems == 1:
            axes = axes.reshape(2, 1)
        
        colors = plt.cm.tab10(np.linspace(0, 1, 10))
        
        for i, (system_name, system_results) in enumerate(results.items()):
            if method not in system_results:
                continue
                
            cluster_data = system_results[method]
            labels = np.array(cluster_data["cluster_labels"])
            
            # Time series of clusters
            time = np.arange(len(labels))  # Frame numbers
            unique_labels = np.unique(labels)
            
            for j, label in enumerate(unique_labels):
                if label == -1:  # Noise in DBSCAN
                    continue
                mask = labels == label
                axes[0, i].scatter(time[mask], np.full(np.sum(mask), label), 
                                 c=[colors[label % len(colors)]], alpha=0.6, 
                                 label=f'Cluster {label}', s=10)
            
            axes[0, i].set_xlabel('Frame')
            axes[0, i].set_ylabel('Cluster ID')
            axes[0, i].set_title(f'{system_name} - Cluster Timeline')
            axes[0, i].legend()
            
            # Cluster population pie chart
            unique_labels_no_noise = unique_labels[unique_labels != -1]
            populations = [np.sum(labels == label) for label in unique_labels_no_noise]
            
            if populations:
                axes[1, i].pie(populations, labels=[f'Cluster {label}' for label in unique_labels_no_noise],
                             autopct='%1.1f%%', colors=[colors[label % len(colors)] for label in unique_labels_no_noise])
                axes[1, i].set_title(f'{system_name} - Cluster Populations')
        
        plt.tight_layout()
        
        if save and self.config.save_plots:
            save_path = self.output_dir / f"clustering_{method}.{self.config.plot_format}"
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Clustering plot saved to {save_path}")
        
        return fig
    
    def plot_contact_map(self, contact_matrix: np.ndarray, system_name: str, save: bool = True) -> plt.Figure:
        """Plot protein contact map.
        
        Args:
            contact_matrix: Contact probability matrix
            system_name: Name of the system
            save: Whether to save the plot
            
        Returns:
            Matplotlib figure
        """
        fig, ax = plt.subplots(figsize=(10, 10))
        
        im = ax.imshow(contact_matrix, cmap='viridis', aspect='equal')
        ax.set_xlabel('Residue Index')
        ax.set_ylabel('Residue Index')
        ax.set_title(f'Contact Map - {system_name}')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Contact Probability')
        
        plt.tight_layout()
        
        if save and self.config.save_plots:
            save_path = self.output_dir / f"contact_map_{system_name}.{self.config.plot_format}"
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Contact map saved to {save_path}")
        
        return fig
    
    def plot_density_profile(self, density_results: Dict[str, Dict], save: bool = True) -> plt.Figure:
        """Plot membrane density profiles.
        
        Args:
            density_results: Dictionary of system_name -> density profile data
            save: Whether to save the plot
            
        Returns:
            Matplotlib figure
        """
        fig, ax = plt.subplots(figsize=(10, 6))
        
        colors = plt.cm.tab10(np.linspace(0, 1, len(density_results)))
        
        for i, (system_name, data) in enumerate(density_results.items()):
            z_centers = data["z_centers"]
            avg_density = data["avg_density"]
            std_density = data["std_density"]
            
            ax.plot(z_centers, avg_density, label=system_name, color=colors[i], linewidth=2)
            ax.fill_between(z_centers, avg_density - std_density, avg_density + std_density,
                           color=colors[i], alpha=0.2)
        
        ax.set_xlabel('Z Position (Å)')
        ax.set_ylabel('Density (atoms/Å³)')
        ax.set_title('Membrane Density Profiles')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save and self.config.save_plots:
            save_path = self.output_dir / f"density_profiles.{self.config.plot_format}"
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Density profiles plot saved to {save_path}")
        
        return fig
    
    def create_interactive_dashboard(self, all_results: Dict[str, Any]) -> go.Figure:
        """Create interactive Plotly dashboard.
        
        Args:
            all_results: All analysis results
            
        Returns:
            Plotly figure
        """
        # Create subplots
        fig = make_subplots(
            rows=3, cols=2,
            subplot_titles=('RMSD', 'Radius of Gyration', 'PIP2 Interactions', 
                          'Membrane Thickness', 'Secondary Structure', 'Clustering'),
            specs=[[{"secondary_y": False}, {"secondary_y": False}],
                   [{"secondary_y": False}, {"secondary_y": False}],
                   [{"secondary_y": False}, {"secondary_y": False}]]
        )
        
        colors = px.colors.qualitative.Set1
        
        # Add traces for each analysis type
        row_col_mapping = {
            'rmsd': (1, 1),
            'radius_of_gyration': (1, 2),
            'pip2_interactions': (2, 1),
            'membrane_thickness': (2, 2),
            'secondary_structure': (3, 1),
            'clustering': (3, 2)
        }
        
        for analysis_type, (row, col) in row_col_mapping.items():
            if analysis_type in all_results:
                for i, (system_name, data) in enumerate(all_results[analysis_type].items()):
                    color = colors[i % len(colors)]
                    
                    if analysis_type == 'rmsd' and isinstance(data, pd.DataFrame):
                        fig.add_trace(
                            go.Scatter(x=data.index/1000, y=data['rmsd_ca'], 
                                     mode='lines', name=f'{system_name} RMSD',
                                     line=dict(color=color)),
                            row=row, col=col
                        )
                    
                    elif analysis_type == 'radius_of_gyration' and isinstance(data, pd.DataFrame):
                        fig.add_trace(
                            go.Scatter(x=data.index/1000, y=data['radius_of_gyration'], 
                                     mode='lines', name=f'{system_name} Rg',
                                     line=dict(color=color)),
                            row=row, col=col
                        )
        
        # Update layout
        fig.update_layout(
            height=900,
            title_text="MD Analysis Dashboard",
            showlegend=True
        )
        
        # Update x-axis labels
        fig.update_xaxes(title_text="Time (ns)", row=3, col=1)
        fig.update_xaxes(title_text="Time (ns)", row=3, col=2)
        
        return fig
    
    def save_interactive_dashboard(self, fig: go.Figure) -> None:
        """Save interactive dashboard as HTML.
        
        Args:
            fig: Plotly figure
        """
        save_path = self.output_dir / "interactive_dashboard.html"
        fig.write_html(str(save_path))
        logger.info(f"Interactive dashboard saved to {save_path}")
    
    def plot_comparison_summary(self, comparison_data: Dict[str, Any], save: bool = True) -> plt.Figure:
        """Plot summary comparison between systems.
        
        Args:
            comparison_data: Comparison results
            save: Whether to save the plot
            
        Returns:
            Matplotlib figure
        """
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        
        # Example comparison plots
        systems = list(comparison_data.keys()) if comparison_data else []
        
        if not systems:
            # Create empty plot
            for ax in axes.flatten():
                ax.text(0.5, 0.5, 'No comparison data available', 
                       ha='center', va='center', transform=ax.transAxes)
                ax.set_xticks([])
                ax.set_yticks([])
        
        plt.suptitle('System Comparison Summary')
        plt.tight_layout()
        
        if save and self.config.save_plots:
            save_path = self.output_dir / f"comparison_summary.{self.config.plot_format}"
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Comparison summary plot saved to {save_path}")
        
        return fig 
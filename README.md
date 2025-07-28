# Membrane Protein Analysis Pipeline

A comprehensive analysis toolkit for molecular dynamics simulations of membrane-protein systems, with a focus on hippocalcin-membrane interactions and PIP2 binding analysis.

## Features

### Core Analysis Capabilities

- **Protein Structure Analysis**
  - RMSD calculations (backbone and C-alpha)
  - Radius of gyration and end-to-end distances
  - Secondary structure evolution (DSSP)
  - Contact map analysis
  - Structural feature extraction for clustering

- **Membrane Analysis**
  - Membrane thickness and area per lipid
  - Lipid composition analysis
  - Membrane curvature and undulation
  - Density profiles along membrane normal
  - Order parameter calculations

- **Protein-Membrane Interactions**
  - Protein-membrane contact analysis
  - Membrane penetration depth
  - PIP2-protein specific interactions
  - Hydrogen bond analysis
  - Electrostatic interaction mapping

- **Clustering Analysis**
  - K-means, DBSCAN, and hierarchical clustering
  - Conformational state identification
  - Cluster transition analysis
  - Comparative clustering between systems

### Advanced Features

- **PIP2 Interaction Analysis**
  - Hydrogen bond detection with configurable cutoffs
  - Electrostatic interaction mapping
  - Binding residue identification
  - Time-resolved interaction dynamics

- **Comparative Analysis**
  - Wild-type vs mutant protein comparison
  - Membrane-only vs protein-containing systems
  - Membrane property perturbation analysis
  - Statistical significance testing

- **Visualization**
  - Matplotlib plots for publication-quality figures
  - Interactive Plotly dashboards
  - Contact maps and density profiles
  - Time series analysis plots

## Installation

### Using Docker (Recommended)

```bash
# Clone the repository
git clone <repository-url>
cd memprot-analysis

# Build Docker image
docker build -t memprot-analysis .

# Run analysis
docker run -v /path/to/your/data:/data memprot-analysis memprot-analyze --help
```

### Using Conda

```bash
# Create conda environment
conda env create -f environment.yml
conda activate memprot

# Install package
pip install -e .
```

### Manual Installation

```bash
# Install dependencies
pip install -r requirements.txt

# Install package
pip install -e .
```

## Quick Start

### 1. Create Configuration File

```bash
memprot-analyze create-config --output my_config.yaml
```

Edit the generated configuration file to specify your simulation files:

```yaml
systems:
  membrane_only:
    topology: "data/membrane_only.gro"
    trajectory: "data/membrane_only.xtc"
    protein_selection: ""
    membrane_selection: "resname POPC POPE POPS CHOL"
    pip2_selection: "resname PIP2"
  
  hippocalcin_wt:
    topology: "data/hippocalcin_wt.gro"
    trajectory: "data/hippocalcin_wt.xtc"
    protein_selection: "protein"
    membrane_selection: "resname POPC POPE POPS CHOL"
    pip2_selection: "resname PIP2"
  
  hippocalcin_n75k:
    topology: "data/hippocalcin_n75k.gro"
    trajectory: "data/hippocalcin_n75k.xtc"
    protein_selection: "protein"
    membrane_selection: "resname POPC POPE POPS CHOL"
    pip2_selection: "resname PIP2"

# Analysis parameters
dt: 10.0  # frame interval in ps
start_time: 1000.0  # equilibration time
end_time: 10000.0  # end time
n_clusters: 5
hbond_cutoff: 3.5
electrostatic_cutoff: 6.0
```

### 2. Run Analysis

```bash
# Full analysis pipeline
memprot-analyze analyze --config my_config.yaml

# Analyze specific systems only
memprot-analyze analyze --config my_config.yaml --systems hippocalcin_wt hippocalcin_n75k

# Skip certain analysis types
memprot-analyze analyze --config my_config.yaml --skip-clustering --skip-plots
```

### 3. Compare Systems

```bash
# Direct comparison between two systems
memprot-analyze compare --config my_config.yaml --system1 hippocalcin_wt --system2 hippocalcin_n75k
```

## Usage Examples

### Python API

```python
from memprot import AnalysisConfig, SimulationSystem
from memprot.analyzers import ProteinAnalyzer, InteractionAnalyzer
from memprot.visualization import PlotManager

# Load configuration
config = AnalysisConfig.from_yaml("config.yaml")

# Load simulation system
system = SimulationSystem(config.systems["hippocalcin_wt"])
system.load()

# Run protein analysis
protein_analyzer = ProteinAnalyzer(config)
protein_analyzer.add_system(system)
protein_analyzer.run_analysis()

# Analyze PIP2 interactions
interaction_analyzer = InteractionAnalyzer(config)
interaction_analyzer.add_system(system)
interaction_analyzer.run_analysis()

# Generate plots
plot_manager = PlotManager(config)
plot_manager.plot_protein_rmsd(protein_analyzer.results['rmsd'])
plot_manager.plot_pip2_interactions(interaction_analyzer.results['pip2_interactions'])
```

### Batch Analysis Script

```python
import yaml
from memprot import AnalysisConfig, SimulationSystem
from memprot.analyzers import ProteinAnalyzer, MembraneAnalyzer, InteractionAnalyzer

# Load configuration
config = AnalysisConfig.from_yaml("batch_config.yaml")

# Load all systems
systems = {}
for name, sys_config in config.systems.items():
    system = SimulationSystem(sys_config)
    system.load()
    systems[name] = system

# Run comparative analysis
analyzers = {}

# Protein analysis
protein_analyzer = ProteinAnalyzer(config)
for system in systems.values():
    protein_analyzer.add_system(system)
protein_analyzer.run_analysis()
analyzers['protein'] = protein_analyzer

# Membrane analysis
membrane_analyzer = MembraneAnalyzer(config)
for system in systems.values():
    membrane_analyzer.add_system(system)
membrane_analyzer.run_analysis()
analyzers['membrane'] = membrane_analyzer

# Calculate membrane influence
membrane_influence = membrane_analyzer.get_membrane_influence_metrics(
    "membrane_only", 
    ["hippocalcin_wt", "hippocalcin_n75k"]
)

print("Membrane Influence Analysis:")
for system, metrics in membrane_influence.items():
    print(f"{system}:")
    for metric, value in metrics.items():
        print(f"  {metric}: {value:.3f}")
```

## Analysis Types

### 1. Protein Structure Analysis

Analyzes protein conformational dynamics:
- RMSD trajectories for stability assessment
- Radius of gyration for compactness
- Secondary structure evolution
- Contact maps for internal interactions

### 2. Membrane Property Analysis

Characterizes membrane properties:
- Thickness variations over time
- Area per lipid calculations
- Curvature and undulation analysis
- Lipid order parameters

### 3. Protein-Membrane Interactions

Studies protein effects on membrane:
- Contact analysis between protein and lipids
- Membrane penetration depth
- Local membrane property perturbations

### 4. PIP2 Interaction Analysis

Focuses on PIP2-protein interactions:
- Hydrogen bond formation and stability
- Electrostatic interactions with basic residues
- Binding site identification
- Interaction dynamics over time

### 5. Clustering Analysis

Identifies conformational states:
- Multiple clustering algorithms (K-means, DBSCAN, hierarchical)
- Optimal cluster number determination
- Transition analysis between states
- Comparative clustering between systems

## Output Files

The analysis generates various output files in the `results/` directory:

```
results/
├── plots/                          # All generated plots
│   ├── protein_rmsd.png
│   ├── pip2_interactions.png
│   ├── membrane_properties.png
│   ├── clustering_kmeans.png
│   └── interactive_dashboard.html
├── proteinanalyzer_results.pkl     # Serialized results
├── interactionanalyzer_results.pkl
├── membraneanalyzer_results.pkl
├── clusteranalyzer_results.pkl
└── analysis_summary.txt            # Summary report
```

## Configuration Options

### System Configuration

```yaml
systems:
  system_name:
    topology: "path/to/topology.gro"
    trajectory: "path/to/trajectory.xtc"
    protein_selection: "protein"  # MDAnalysis selection
    membrane_selection: "resname POPC POPE POPS CHOL"
    pip2_selection: "resname PIP2"
```

### Analysis Parameters

```yaml
# Time settings
dt: 10.0                # Frame interval (ps)
start_time: 1000.0      # Analysis start time (ps)
end_time: 10000.0       # Analysis end time (ps)

# Clustering settings
clustering_method: "kmeans"  # or "dbscan", "hierarchical", "all"
n_clusters: 5
clustering_features:
  - "radius_of_gyration"
  - "end_to_end_distance"
  - "secondary_structure"

# Interaction analysis
hbond_cutoff: 3.5           # Hydrogen bond distance cutoff (Å)
hbond_angle_cutoff: 150.0   # Hydrogen bond angle cutoff (degrees)
electrostatic_cutoff: 6.0   # Electrostatic interaction cutoff (Å)

# Output settings
output_dir: "results"
save_plots: true
plot_format: "png"         # or "pdf", "svg"
n_jobs: -1                 # Parallel processing cores
```

## Typical Analysis Workflow

### For Hippocalcin-Membrane Systems

1. **Prepare three systems:**
   - Membrane only (baseline)
   - Wild-type hippocalcin + membrane
   - N75K mutant hippocalcin + membrane

2. **Run comprehensive analysis:**
   ```bash
   memprot-analyze analyze --config hippocalcin_config.yaml
   ```

3. **Key questions answered:**
   - How does protein binding affect membrane properties?
   - What are the differences between wild-type and mutant?
   - Which residues are critical for PIP2 binding?
   - How stable are the protein conformations?

4. **Examine results:**
   - Check `analysis_summary.txt` for overview
   - View plots in `results/plots/`
   - Use interactive dashboard for detailed exploration

### Specific Analysis Examples

#### PIP2 Binding Analysis
```python
# Get PIP2 binding residues
binding_residues = interaction_analyzer.get_pip2_binding_residues(
    "hippocalcin_wt", 
    contact_threshold=0.3
)
print(f"PIP2 binding residues: {binding_residues}")

# Compare binding between systems
comparison = interaction_analyzer.compare_pip2_interactions(
    "hippocalcin_wt", 
    "hippocalcin_n75k"
)
```

#### Membrane Influence Assessment
```python
# Calculate membrane perturbation
influence_metrics = membrane_analyzer.get_membrane_influence_metrics(
    "membrane_only",
    ["hippocalcin_wt", "hippocalcin_n75k"]
)

# Print thickness changes
for system, metrics in influence_metrics.items():
    print(f"{system} thickness change: {metrics['thickness_change']:.2f} Å")
```

#### Conformational Clustering
```python
# Get cluster transitions
transitions = cluster_analyzer.get_cluster_transitions("hippocalcin_wt")
print("Cluster transition matrix:")
print(transitions)

# Compare clustering between systems
comparison = cluster_analyzer.compare_systems_clustering(
    "hippocalcin_wt", 
    "hippocalcin_n75k"
)
```

## Dependencies

### Core Dependencies
- MDAnalysis >= 2.4.0
- NumPy >= 1.21.0
- Pandas >= 1.3.0
- SciPy >= 1.7.0
- Scikit-learn >= 1.0.0

### Visualization
- Matplotlib >= 3.4.0
- Seaborn >= 0.11.0
- Plotly >= 5.0.0

### Optional
- Jupyter for interactive analysis
- NGLView for 3D structure visualization

## Contributing

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## License

MIT License - see LICENSE file for details.

## Citation

If you use this software in your research, please cite:

```
@software{memprot_analysis,
  title={Membrane Protein Analysis Pipeline},
  author={MD Analysis Team},
  year={2024},
  url={https://github.com/example/memprot-analysis}
}
```

## Support

For questions and support:
- Check the documentation
- Open an issue on GitHub
- Contact the development team

## Changelog

### Version 1.0.0
- Initial release
- Complete analysis pipeline for hippocalcin-membrane systems
- PIP2 interaction analysis
- Clustering and comparative analysis
- Docker containerization
- Interactive visualization dashboard 
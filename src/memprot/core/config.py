"""Configuration management for MD analysis pipeline."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Union
import yaml
from loguru import logger


@dataclass
class SystemConfig:
    """Configuration for simulation system files."""
    topology: Path
    trajectory: Path
    name: str
    protein_selection: str = "protein"
    membrane_selection: str = "resname POPC POPE POPS CHOL"
    pip2_selection: str = "resname PIP2"
    
    def __post_init__(self):
        self.topology = Path(self.topology)
        self.trajectory = Path(self.trajectory)


@dataclass
class AnalysisConfig:
    """Main configuration class for MD analysis."""
    
    # System configurations
    systems: Dict[str, SystemConfig] = field(default_factory=dict)
    
    # Analysis parameters
    dt: float = 1.0  # ps, frame interval
    start_time: float = 0.0  # ps
    end_time: Optional[float] = None  # ps
    
    # Clustering parameters
    clustering_method: str = "kmeans"
    n_clusters: int = 5
    clustering_features: List[str] = field(default_factory=lambda: [
        "radius_of_gyration", "end_to_end_distance", "secondary_structure"
    ])
    
    # Interaction analysis
    hbond_cutoff: float = 3.5  # Angstrom
    hbond_angle_cutoff: float = 150.0  # degrees
    electrostatic_cutoff: float = 6.0  # Angstrom
    
    # Output settings
    output_dir: Path = field(default_factory=lambda: Path("results"))
    save_plots: bool = True
    plot_format: str = "png"
    
    # Parallel processing
    n_jobs: int = -1
    
    @classmethod
    def from_yaml(cls, config_file: Union[str, Path]) -> "AnalysisConfig":
        """Load configuration from YAML file."""
        config_file = Path(config_file)
        
        if not config_file.exists():
            raise FileNotFoundError(f"Config file not found: {config_file}")
            
        with open(config_file, "r") as f:
            data = yaml.safe_load(f)
            
        # Convert system configurations
        systems = {}
        if "systems" in data:
            for name, sys_data in data["systems"].items():
                systems[name] = SystemConfig(name=name, **sys_data)
        data["systems"] = systems
        
        # Convert paths
        if "output_dir" in data:
            data["output_dir"] = Path(data["output_dir"])
            
        return cls(**data)
    
    def to_yaml(self, config_file: Union[str, Path]) -> None:
        """Save configuration to YAML file."""
        config_file = Path(config_file)
        config_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Convert to serializable format
        data = {}
        for field_name, field_value in self.__dict__.items():
            if field_name == "systems":
                data[field_name] = {
                    name: {
                        "topology": str(sys.topology),
                        "trajectory": str(sys.trajectory),
                        "protein_selection": sys.protein_selection,
                        "membrane_selection": sys.membrane_selection,
                        "pip2_selection": sys.pip2_selection,
                    }
                    for name, sys in field_value.items()
                }
            elif isinstance(field_value, Path):
                data[field_name] = str(field_value)
            else:
                data[field_name] = field_value
                
        with open(config_file, "w") as f:
            yaml.dump(data, f, default_flow_style=False, indent=2)
            
        logger.info(f"Configuration saved to {config_file}")
    
    def validate(self) -> None:
        """Validate configuration parameters."""
        errors = []
        
        # Check systems
        if not self.systems:
            errors.append("No systems defined in configuration")
            
        for name, system in self.systems.items():
            if not system.topology.exists():
                errors.append(f"Topology file not found for {name}: {system.topology}")
            if not system.trajectory.exists():
                errors.append(f"Trajectory file not found for {name}: {system.trajectory}")
                
        # Check analysis parameters
        if self.dt <= 0:
            errors.append("dt must be positive")
            
        if self.n_clusters <= 0:
            errors.append("n_clusters must be positive")
            
        if self.hbond_cutoff <= 0:
            errors.append("hbond_cutoff must be positive")
            
        if errors:
            raise ValueError("Configuration validation failed:\n" + "\n".join(errors))
            
        logger.info("Configuration validation passed")


def create_example_config() -> AnalysisConfig:
    """Create an example configuration."""
    config = AnalysisConfig()
    
    # Add example systems
    config.systems = {
        "membrane_only": SystemConfig(
            name="membrane_only",
            topology=Path("data/membrane_only.gro"),
            trajectory=Path("data/membrane_only.xtc"),
            protein_selection="",
            membrane_selection="resname POPC POPE POPS CHOL",
            pip2_selection="resname PIP2"
        ),
        "hippocalcin_wt": SystemConfig(
            name="hippocalcin_wt",
            topology=Path("data/hippocalcin_wt.gro"), 
            trajectory=Path("data/hippocalcin_wt.xtc"),
            protein_selection="protein",
            membrane_selection="resname POPC POPE POPS CHOL",
            pip2_selection="resname PIP2"
        ),
        "hippocalcin_n75k": SystemConfig(
            name="hippocalcin_n75k",
            topology=Path("data/hippocalcin_n75k.gro"),
            trajectory=Path("data/hippocalcin_n75k.xtc"), 
            protein_selection="protein",
            membrane_selection="resname POPC POPE POPS CHOL",
            pip2_selection="resname PIP2"
        )
    }
    
    return config 
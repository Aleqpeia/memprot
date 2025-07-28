"""Base analyzer class with common functionality."""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
import numpy as np
import pandas as pd
from loguru import logger

from .system import SimulationSystem
from .config import AnalysisConfig


class BaseAnalyzer(ABC):
    """Base class for all analyzers."""
    
    def __init__(self, config: AnalysisConfig):
        """Initialize analyzer.
        
        Args:
            config: Analysis configuration
        """
        self.config = config
        self.results: Dict[str, Any] = {}
        self._systems: Dict[str, SimulationSystem] = {}
        
    def add_system(self, system: SimulationSystem) -> None:
        """Add a simulation system.
        
        Args:
            system: Simulation system to add
        """
        self._systems[system.config.name] = system
        logger.info(f"Added system: {system.config.name}")
        
    def get_system(self, name: str) -> SimulationSystem:
        """Get simulation system by name.
        
        Args:
            name: System name
            
        Returns:
            Simulation system
        """
        if name not in self._systems:
            raise KeyError(f"System '{name}' not found")
        return self._systems[name]
    
    def get_frame_slice(self, system_name: str) -> range:
        """Get frame slice for analysis.
        
        Args:
            system_name: Name of system
            
        Returns:
            Range of frame indices
        """
        system = self.get_system(system_name)
        return system.get_trajectory_slice(
            start_time=self.config.start_time,
            end_time=self.config.end_time,
            dt=self.config.dt
        )
    
    @abstractmethod
    def run_analysis(self, system_names: Optional[List[str]] = None) -> None:
        """Run analysis on specified systems.
        
        Args:
            system_names: Names of systems to analyze. If None, analyze all.
        """
        pass
    
    def save_results(self, output_dir: Optional[Path] = None) -> None:
        """Save analysis results.
        
        Args:
            output_dir: Output directory. Uses config default if None.
        """
        if output_dir is None:
            output_dir = self.config.output_dir
            
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save results as pickle and JSON where possible
        results_file = output_dir / f"{self.__class__.__name__.lower()}_results.pkl"
        
        try:
            import pickle
            with open(results_file, "wb") as f:
                pickle.dump(self.results, f)
            logger.info(f"Results saved to {results_file}")
        except Exception as e:
            logger.error(f"Failed to save results: {e}")
    
    def _validate_systems(self, system_names: List[str]) -> None:
        """Validate that required systems are available.
        
        Args:
            system_names: Names of systems to validate
        """
        missing_systems = [name for name in system_names if name not in self._systems]
        if missing_systems:
            raise ValueError(f"Missing systems: {missing_systems}")
    
    def _get_time_array(self, system_name: str) -> np.ndarray:
        """Get time array for frames.
        
        Args:
            system_name: Name of system
            
        Returns:
            Array of times in ps
        """
        system = self.get_system(system_name)
        frame_slice = self.get_frame_slice(system_name)
        
        times = []
        for frame_idx in frame_slice:
            system.universe.trajectory[frame_idx]
            times.append(system.universe.trajectory.time)
            
        return np.array(times)
    
    def _create_dataframe(
        self, 
        data: Dict[str, np.ndarray], 
        system_name: str
    ) -> pd.DataFrame:
        """Create DataFrame with time index.
        
        Args:
            data: Dictionary of arrays with analysis data
            system_name: Name of system
            
        Returns:
            DataFrame with time as index
        """
        times = self._get_time_array(system_name)
        
        # Ensure all arrays have same length as times
        for key, array in data.items():
            if len(array) != len(times):
                raise ValueError(
                    f"Data array '{key}' length {len(array)} != times length {len(times)}"
                )
        
        df = pd.DataFrame(data)
        df.index = times
        df.index.name = "time_ps"
        
        return df
    
    def get_summary_statistics(self) -> Dict[str, Any]:
        """Get summary statistics for all results.
        
        Returns:
            Dictionary of summary statistics
        """
        summary = {}
        
        for key, data in self.results.items():
            if isinstance(data, (pd.DataFrame, pd.Series)):
                if isinstance(data, pd.DataFrame):
                    summary[key] = {
                        col: {
                            "mean": float(data[col].mean()),
                            "std": float(data[col].std()),
                            "min": float(data[col].min()),
                            "max": float(data[col].max()),
                            "median": float(data[col].median()),
                        }
                        for col in data.select_dtypes(include=[np.number]).columns
                    }
                else:
                    if data.dtype in [np.float64, np.float32, np.int64, np.int32]:
                        summary[key] = {
                            "mean": float(data.mean()),
                            "std": float(data.std()),
                            "min": float(data.min()),
                            "max": float(data.max()),
                            "median": float(data.median()),
                        }
            elif isinstance(data, np.ndarray) and data.dtype in [np.float64, np.float32]:
                summary[key] = {
                    "mean": float(np.mean(data)),
                    "std": float(np.std(data)),
                    "min": float(np.min(data)),
                    "max": float(np.max(data)),
                    "median": float(np.median(data)),
                }
                
        return summary 
"""Simulation system management and data loading."""

from pathlib import Path
from typing import Optional, Union, Dict, Any
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from loguru import logger

from .config import SystemConfig


class SimulationSystem:
    """Manages loading and basic operations for MD simulation data."""
    
    def __init__(self, config: SystemConfig):
        """Initialize simulation system.
        
        Args:
            config: System configuration containing file paths and selections
        """
        self.config = config
        self.universe: Optional[mda.Universe] = None
        self._protein: Optional[mda.AtomGroup] = None
        self._membrane: Optional[mda.AtomGroup] = None
        self._pip2: Optional[mda.AtomGroup] = None
        self._loaded = False
        
    def load(self) -> None:
        """Load topology and trajectory files."""
        try:
            logger.info(f"Loading system: {self.config.name}")
            logger.info(f"Topology: {self.config.topology}")
            logger.info(f"Trajectory: {self.config.trajectory}")
            
            self.universe = mda.Universe(
                str(self.config.topology),
                str(self.config.trajectory)
            )
            
            logger.info(f"Loaded {len(self.universe.trajectory)} frames")
            logger.info(f"System contains {len(self.universe.atoms)} atoms")
            
            # Create atom selections
            self._create_selections()
            self._loaded = True
            
        except Exception as e:
            logger.error(f"Failed to load system {self.config.name}: {e}")
            raise
    
    def _create_selections(self) -> None:
        """Create atom group selections."""
        if not self.universe:
            raise RuntimeError("Universe not loaded")
            
        # Protein selection
        if self.config.protein_selection:
            try:
                self._protein = self.universe.select_atoms(self.config.protein_selection)
                logger.info(f"Protein selection: {len(self._protein)} atoms")
            except Exception as e:
                logger.warning(f"Failed to create protein selection: {e}")
                
        # Membrane selection  
        if self.config.membrane_selection:
            try:
                self._membrane = self.universe.select_atoms(self.config.membrane_selection)
                logger.info(f"Membrane selection: {len(self._membrane)} atoms")
            except Exception as e:
                logger.warning(f"Failed to create membrane selection: {e}")
                
        # PIP2 selection
        if self.config.pip2_selection:
            try:
                self._pip2 = self.universe.select_atoms(self.config.pip2_selection)
                logger.info(f"PIP2 selection: {len(self._pip2)} atoms")
            except Exception as e:
                logger.warning(f"Failed to create PIP2 selection: {e}")
    
    @property
    def protein(self) -> Optional[mda.AtomGroup]:
        """Get protein atom group."""
        if not self._loaded:
            raise RuntimeError("System not loaded")
        return self._protein
    
    @property
    def membrane(self) -> Optional[mda.AtomGroup]:
        """Get membrane atom group."""
        if not self._loaded:
            raise RuntimeError("System not loaded")
        return self._membrane
    
    @property 
    def pip2(self) -> Optional[mda.AtomGroup]:
        """Get PIP2 atom group."""
        if not self._loaded:
            raise RuntimeError("System not loaded")
        return self._pip2
    
    @property
    def has_protein(self) -> bool:
        """Check if system has protein."""
        return self._protein is not None and len(self._protein) > 0
    
    @property
    def has_membrane(self) -> bool:
        """Check if system has membrane."""
        return self._membrane is not None and len(self._membrane) > 0
        
    @property
    def has_pip2(self) -> bool:
        """Check if system has PIP2."""
        return self._pip2 is not None and len(self._pip2) > 0
    
    def get_trajectory_slice(
        self, 
        start_time: Optional[float] = None,
        end_time: Optional[float] = None,
        dt: Optional[float] = None
    ) -> range:
        """Get frame indices for trajectory slice.
        
        Args:
            start_time: Start time in ps
            end_time: End time in ps  
            dt: Frame interval in ps
            
        Returns:
            Range of frame indices
        """
        if not self.universe:
            raise RuntimeError("Universe not loaded")
            
        traj = self.universe.trajectory
        total_time = traj.totaltime
        frame_dt = traj.dt
        
        start_frame = 0
        if start_time is not None:
            start_frame = max(0, int(start_time / frame_dt))
            
        end_frame = len(traj)
        if end_time is not None:
            end_frame = min(len(traj), int(end_time / frame_dt))
            
        step = 1
        if dt is not None and dt > frame_dt:
            step = int(dt / frame_dt)
            
        return range(start_frame, end_frame, step)
    
    def align_protein_trajectory(
        self, 
        reference_frame: int = 0,
        selection: str = "name CA"
    ) -> None:
        """Align protein trajectory to reference frame.
        
        Args:
            reference_frame: Frame to use as reference
            selection: Atoms to use for alignment
        """
        if not self.has_protein:
            logger.warning("No protein found, skipping alignment")
            return
            
        logger.info(f"Aligning protein trajectory to frame {reference_frame}")
        
        try:
            # Create alignment selection
            align_atoms = self.protein.select_atoms(selection)
            
            # Set reference
            self.universe.trajectory[reference_frame]
            reference_coords = align_atoms.positions.copy()
            
            # Align trajectory
            align.AlignTraj(
                self.universe,
                self.universe,
                select=f"({self.config.protein_selection}) and ({selection})",
                filename=None,
                in_memory=True
            ).run()
            
            logger.info("Protein trajectory aligned successfully")
            
        except Exception as e:
            logger.error(f"Failed to align trajectory: {e}")
            raise
    
    def get_system_info(self) -> Dict[str, Any]:
        """Get system information summary."""
        if not self._loaded:
            return {"loaded": False}
            
        info = {
            "loaded": True,
            "name": self.config.name,
            "n_frames": len(self.universe.trajectory),
            "n_atoms": len(self.universe.atoms),
            "total_time_ps": self.universe.trajectory.totaltime,
            "frame_dt_ps": self.universe.trajectory.dt,
            "has_protein": self.has_protein,
            "has_membrane": self.has_membrane, 
            "has_pip2": self.has_pip2,
        }
        
        if self.has_protein:
            info["n_protein_atoms"] = len(self._protein)
            info["n_residues"] = len(self._protein.residues)
            
        if self.has_membrane:
            info["n_membrane_atoms"] = len(self._membrane)
            
        if self.has_pip2:
            info["n_pip2_atoms"] = len(self._pip2)
            info["n_pip2_molecules"] = len(self._pip2.residues)
            
        return info
    
    def __repr__(self) -> str:
        """String representation."""
        if not self._loaded:
            return f"SimulationSystem('{self.config.name}', not loaded)"
            
        return (
            f"SimulationSystem('{self.config.name}', "
            f"{len(self.universe.trajectory)} frames, "
            f"{len(self.universe.atoms)} atoms)"
        ) 
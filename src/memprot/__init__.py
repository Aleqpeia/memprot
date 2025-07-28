"""
Membrane Protein Analysis Package
=================================

A comprehensive toolkit for analyzing molecular dynamics simulations
of membrane-protein systems, with focus on hippocalcin interactions.
"""

__version__ = "1.0.0"
__author__ = "MD Analysis Team"

from .core import SimulationSystem, AnalysisConfig
from .analyzers import (
    ProteinAnalyzer,
    MembraneAnalyzer, 
    InteractionAnalyzer,
    ClusterAnalyzer,
)
from .visualization import PlotManager
from .utils import setup_logging

__all__ = [
    "SimulationSystem",
    "AnalysisConfig", 
    "ProteinAnalyzer",
    "MembraneAnalyzer",
    "InteractionAnalyzer", 
    "ClusterAnalyzer",
    "PlotManager",
    "setup_logging",
] 
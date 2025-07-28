"""Clustering analysis for protein conformations."""

from typing import Dict, List, Optional, Tuple, Any
import numpy as np
import pandas as pd
from loguru import logger
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score, calinski_harabasz_score

from ..core.base import BaseAnalyzer


class ClusterAnalyzer(BaseAnalyzer):
    """Analyzes protein conformational clusters."""
    
    def run_analysis(self, system_names: Optional[List[str]] = None) -> None:
        """Run clustering analysis.
        
        Args:
            system_names: Names of systems to analyze. If None, analyze all with protein.
        """
        if system_names is None:
            system_names = [name for name, sys in self._systems.items() if sys.has_protein]
        else:
            system_names = [name for name in system_names 
                           if name in self._systems and self._systems[name].has_protein]
        
        if not system_names:
            logger.warning("No systems with protein found")
            return
            
        self._validate_systems(system_names)
        
        # Get structural features from ProteinAnalyzer
        self._prepare_clustering_data(system_names)
        
        # Perform clustering for each system
        for system_name in system_names:
            logger.info(f"Clustering analysis for system: {system_name}")
            self._cluster_system(system_name)
        
        # Compare clusters between systems
        if len(system_names) > 1:
            self._compare_clusters(system_names)
            
        logger.info("Clustering analysis completed")
    
    def _prepare_clustering_data(self, system_names: List[str]) -> None:
        """Prepare data for clustering analysis."""
        if "clustering_data" not in self.results:
            self.results["clustering_data"] = {}
            
        for system_name in system_names:
            system = self.get_system(system_name)
            
            # Collect features based on configuration
            features = {}
            
            # Try to get structural features from ProteinAnalyzer results
            # This assumes ProteinAnalyzer has been run first
            protein_results = getattr(self, '_protein_analyzer_results', None)
            
            if protein_results is None:
                # Generate basic features directly
                features = self._generate_basic_features(system_name)
            else:
                # Use results from ProteinAnalyzer
                features = self._extract_protein_features(system_name, protein_results)
            
            if features:
                self.results["clustering_data"][system_name] = features
                logger.info(f"Prepared {len(features)} features for {system_name}")
            else:
                logger.warning(f"No clustering features available for {system_name}")
    
    def _generate_basic_features(self, system_name: str) -> Dict[str, np.ndarray]:
        """Generate basic structural features for clustering."""
        system = self.get_system(system_name)
        protein = system.protein
        
        if protein is None:
            return {}
            
        frame_slice = self.get_frame_slice(system_name)
        
        features = {}
        
        # Radius of gyration
        if "radius_of_gyration" in self.config.clustering_features:
            rg_values = []
            for frame_idx in frame_slice:
                system.universe.trajectory[frame_idx]
                rg_values.append(protein.radius_of_gyration())
            features["radius_of_gyration"] = np.array(rg_values)
        
        # End-to-end distance
        if "end_to_end_distance" in self.config.clustering_features:
            residues = protein.residues
            n_terminus = residues[0].atoms.select_atoms("name CA")[0]
            c_terminus = residues[-1].atoms.select_atoms("name CA")[0]
            
            distances = []
            for frame_idx in frame_slice:
                system.universe.trajectory[frame_idx]
                distance = np.linalg.norm(c_terminus.position - n_terminus.position)
                distances.append(distance)
            features["end_to_end_distance"] = np.array(distances)
        
        # Principal components of atomic coordinates
        if "pca_coords" in self.config.clustering_features:
            ca_atoms = protein.select_atoms("name CA")
            coords_array = []
            
            for frame_idx in frame_slice:
                system.universe.trajectory[frame_idx]
                coords_array.append(ca_atoms.positions.flatten())
            
            coords_array = np.array(coords_array)
            
            # Perform PCA
            pca = PCA(n_components=min(10, coords_array.shape[1]))
            pca_coords = pca.fit_transform(coords_array)
            
            for i in range(pca_coords.shape[1]):
                features[f"pca_{i+1}"] = pca_coords[:, i]
        
        return features
    
    def _extract_protein_features(self, system_name: str, protein_results: Dict) -> Dict[str, np.ndarray]:
        """Extract features from ProteinAnalyzer results."""
        features = {}
        
        # Radius of gyration
        if ("radius_of_gyration" in protein_results and 
            system_name in protein_results["radius_of_gyration"]):
            rg_data = protein_results["radius_of_gyration"][system_name]
            features["radius_of_gyration"] = rg_data["radius_of_gyration"].values
        
        # End-to-end distance
        if ("end_to_end_distance" in protein_results and 
            system_name in protein_results["end_to_end_distance"]):
            ete_data = protein_results["end_to_end_distance"][system_name]
            features["end_to_end_distance"] = ete_data["end_to_end_distance"].values
        
        # Secondary structure
        if ("secondary_structure" in protein_results and 
            system_name in protein_results["secondary_structure"]):
            ss_data = protein_results["secondary_structure"][system_name]
            features["helix_frac"] = ss_data["helix_frac"].values
            features["sheet_frac"] = ss_data["sheet_frac"].values
            features["coil_frac"] = ss_data["coil_frac"].values
        
        return features
    
    def _cluster_system(self, system_name: str) -> None:
        """Perform clustering for a single system."""
        if system_name not in self.results["clustering_data"]:
            logger.warning(f"No clustering data for {system_name}")
            return
            
        features_dict = self.results["clustering_data"][system_name]
        
        # Combine features into matrix
        feature_names = list(features_dict.keys())
        feature_matrix = np.column_stack([features_dict[name] for name in feature_names])
        
        # Standardize features
        scaler = StandardScaler()
        feature_matrix_scaled = scaler.fit_transform(feature_matrix)
        
        # Store preprocessing info
        if "clustering_preprocessing" not in self.results:
            self.results["clustering_preprocessing"] = {}
        self.results["clustering_preprocessing"][system_name] = {
            "feature_names": feature_names,
            "scaler": scaler
        }
        
        # Perform different clustering methods
        clustering_results = {}
        
        # K-means clustering
        if self.config.clustering_method in ["kmeans", "all"]:
            kmeans_results = self._perform_kmeans(feature_matrix_scaled, system_name)
            clustering_results["kmeans"] = kmeans_results
        
        # DBSCAN clustering
        if self.config.clustering_method in ["dbscan", "all"]:
            dbscan_results = self._perform_dbscan(feature_matrix_scaled, system_name)
            clustering_results["dbscan"] = dbscan_results
        
        # Hierarchical clustering
        if self.config.clustering_method in ["hierarchical", "all"]:
            hierarchical_results = self._perform_hierarchical(feature_matrix_scaled, system_name)
            clustering_results["hierarchical"] = hierarchical_results
        
        # Store results
        if "clustering_results" not in self.results:
            self.results["clustering_results"] = {}
        self.results["clustering_results"][system_name] = clustering_results
        
        # Analyze clusters
        self._analyze_clusters(system_name, feature_matrix_scaled, clustering_results)
    
    def _perform_kmeans(self, feature_matrix: np.ndarray, system_name: str) -> Dict[str, Any]:
        """Perform K-means clustering."""
        results = {}
        
        n_clusters = self.config.n_clusters
        
        # Try different numbers of clusters
        cluster_range = range(2, min(n_clusters + 3, feature_matrix.shape[0] // 2))
        
        silhouette_scores = []
        calinski_scores = []
        inertias = []
        
        for n in cluster_range:
            kmeans = KMeans(n_clusters=n, random_state=42, n_init=10)
            labels = kmeans.fit_predict(feature_matrix)
            
            if len(np.unique(labels)) > 1:
                sil_score = silhouette_score(feature_matrix, labels)
                cal_score = calinski_harabasz_score(feature_matrix, labels)
            else:
                sil_score = -1
                cal_score = 0
            
            silhouette_scores.append(sil_score)
            calinski_scores.append(cal_score)
            inertias.append(kmeans.inertia_)
        
        # Choose best number of clusters based on silhouette score
        best_n = cluster_range[np.argmax(silhouette_scores)]
        
        # Final clustering with best n
        final_kmeans = KMeans(n_clusters=best_n, random_state=42, n_init=10)
        final_labels = final_kmeans.fit_predict(feature_matrix)
        
        results = {
            "labels": final_labels,
            "n_clusters": best_n,
            "centers": final_kmeans.cluster_centers_,
            "silhouette_score": max(silhouette_scores),
            "calinski_score": calinski_scores[np.argmax(silhouette_scores)],
            "inertia": final_kmeans.inertia_,
            "cluster_evaluation": {
                "n_clusters_range": list(cluster_range),
                "silhouette_scores": silhouette_scores,
                "calinski_scores": calinski_scores,
                "inertias": inertias
            }
        }
        
        return results
    
    def _perform_dbscan(self, feature_matrix: np.ndarray, system_name: str) -> Dict[str, Any]:
        """Perform DBSCAN clustering."""
        # Try different epsilon values
        eps_values = np.linspace(0.1, 2.0, 10)
        min_samples = max(2, feature_matrix.shape[0] // 50)
        
        best_eps = None
        best_score = -1
        best_labels = None
        
        for eps in eps_values:
            dbscan = DBSCAN(eps=eps, min_samples=min_samples)
            labels = dbscan.fit_predict(feature_matrix)
            
            n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
            
            if n_clusters > 1:
                try:
                    score = silhouette_score(feature_matrix, labels)
                    if score > best_score:
                        best_score = score
                        best_eps = eps
                        best_labels = labels
                except:
                    continue
        
        if best_labels is None:
            # Fallback to simple DBSCAN
            dbscan = DBSCAN(eps=0.5, min_samples=min_samples)
            best_labels = dbscan.fit_predict(feature_matrix)
            best_eps = 0.5
            best_score = -1
        
        results = {
            "labels": best_labels,
            "n_clusters": len(set(best_labels)) - (1 if -1 in best_labels else 0),
            "eps": best_eps,
            "min_samples": min_samples,
            "silhouette_score": best_score,
            "n_noise": np.sum(best_labels == -1)
        }
        
        return results
    
    def _perform_hierarchical(self, feature_matrix: np.ndarray, system_name: str) -> Dict[str, Any]:
        """Perform hierarchical clustering."""
        n_clusters = self.config.n_clusters
        
        clustering = AgglomerativeClustering(
            n_clusters=n_clusters,
            linkage='ward'
        )
        labels = clustering.fit_predict(feature_matrix)
        
        silhouette = silhouette_score(feature_matrix, labels)
        calinski = calinski_harabasz_score(feature_matrix, labels)
        
        results = {
            "labels": labels,
            "n_clusters": n_clusters,
            "silhouette_score": silhouette,
            "calinski_score": calinski
        }
        
        return results
    
    def _analyze_clusters(self, system_name: str, feature_matrix: np.ndarray, clustering_results: Dict) -> None:
        """Analyze cluster properties."""
        if "cluster_analysis" not in self.results:
            self.results["cluster_analysis"] = {}
        
        analysis_results = {}
        
        for method, results in clustering_results.items():
            labels = results["labels"]
            unique_labels = np.unique(labels)
            
            cluster_stats = {}
            
            for cluster_id in unique_labels:
                if cluster_id == -1:  # Noise in DBSCAN
                    continue
                    
                cluster_mask = labels == cluster_id
                cluster_features = feature_matrix[cluster_mask]
                
                stats = {
                    "size": int(np.sum(cluster_mask)),
                    "fraction": float(np.mean(cluster_mask)),
                    "center": np.mean(cluster_features, axis=0).tolist(),
                    "std": np.std(cluster_features, axis=0).tolist()
                }
                
                cluster_stats[f"cluster_{cluster_id}"] = stats
            
            analysis_results[method] = {
                "cluster_stats": cluster_stats,
                "cluster_labels": labels.tolist()
            }
        
        self.results["cluster_analysis"][system_name] = analysis_results
    
    def _compare_clusters(self, system_names: List[str]) -> None:
        """Compare clusters between different systems."""
        if "cluster_comparison" not in self.results:
            self.results["cluster_comparison"] = {}
        
        # Compare cluster populations
        comparison_results = {}
        
        for method in ["kmeans", "dbscan", "hierarchical"]:
            method_comparison = {}
            
            for system_name in system_names:
                if (system_name in self.results.get("cluster_analysis", {}) and
                    method in self.results["cluster_analysis"][system_name]):
                    
                    cluster_stats = self.results["cluster_analysis"][system_name][method]["cluster_stats"]
                    
                    # Extract cluster populations
                    populations = {}
                    for cluster_name, stats in cluster_stats.items():
                        populations[cluster_name] = stats["size"]
                    
                    method_comparison[system_name] = populations
            
            if len(method_comparison) > 1:
                comparison_results[method] = method_comparison
        
        self.results["cluster_comparison"] = comparison_results
    
    def get_cluster_transitions(self, system_name: str, method: str = "kmeans") -> pd.DataFrame:
        """Get cluster transition matrix.
        
        Args:
            system_name: Name of system
            method: Clustering method
            
        Returns:
            Transition probability matrix
        """
        if (system_name not in self.results.get("cluster_analysis", {}) or
            method not in self.results["cluster_analysis"][system_name]):
            raise ValueError(f"No {method} clustering results for {system_name}")
        
        labels = self.results["cluster_analysis"][system_name][method]["cluster_labels"]
        labels = np.array(labels)
        
        unique_labels = np.unique(labels)
        n_clusters = len(unique_labels)
        
        # Create transition matrix
        transitions = np.zeros((n_clusters, n_clusters))
        
        for i in range(len(labels) - 1):
            from_cluster = np.where(unique_labels == labels[i])[0][0]
            to_cluster = np.where(unique_labels == labels[i + 1])[0][0]
            transitions[from_cluster, to_cluster] += 1
        
        # Normalize to probabilities
        row_sums = transitions.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1  # Avoid division by zero
        transitions = transitions / row_sums
        
        # Create DataFrame
        df = pd.DataFrame(
            transitions,
            index=[f"cluster_{label}" for label in unique_labels],
            columns=[f"cluster_{label}" for label in unique_labels]
        )
        
        return df
    
    def compare_systems_clustering(self, system1: str, system2: str, method: str = "kmeans") -> Dict[str, Any]:
        """Compare clustering results between two systems.
        
        Args:
            system1: First system name
            system2: Second system name
            method: Clustering method
            
        Returns:
            Comparison results
        """
        comparison = {}
        
        # Get cluster statistics
        if (system1 in self.results.get("cluster_analysis", {}) and
            system2 in self.results.get("cluster_analysis", {}) and
            method in self.results["cluster_analysis"][system1] and
            method in self.results["cluster_analysis"][system2]):
            
            stats1 = self.results["cluster_analysis"][system1][method]["cluster_stats"]
            stats2 = self.results["cluster_analysis"][system2][method]["cluster_stats"]
            
            # Compare cluster populations
            pop_comparison = {}
            all_clusters = set(stats1.keys()) | set(stats2.keys())
            
            for cluster in all_clusters:
                size1 = stats1.get(cluster, {}).get("size", 0)
                size2 = stats2.get(cluster, {}).get("size", 0)
                frac1 = stats1.get(cluster, {}).get("fraction", 0)
                frac2 = stats2.get(cluster, {}).get("fraction", 0)
                
                pop_comparison[cluster] = {
                    f"{system1}_size": size1,
                    f"{system2}_size": size2,
                    f"{system1}_fraction": frac1,
                    f"{system2}_fraction": frac2,
                    "size_ratio": size2 / max(size1, 1),
                    "fraction_diff": frac2 - frac1
                }
            
            comparison["population_comparison"] = pop_comparison
        
        return comparison 
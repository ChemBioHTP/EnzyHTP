"""This submodule contains code for cluster representative structures from a structure ensemble
+ rep_stru_clustering()
    Cluster the structure ensemble into representative structure

Author: Qianzhen Shao <shaoqz@icloud.com>

Date: 2025-07-08
"""
from typing import Callable, Dict, Tuple
import numpy as np
from sklearn.cluster import KMeans

from enzy_htp.preparation import remove_solvent
from enzy_htp import interface as eh_interface
from enzy_htp.structure import StruSelection, Structure
from enzy_htp.structure.structure_ensemble import StructureEnsemble
from enzy_htp.structure.structure_selection import select_stru

from .rmsd import rmsd

def rep_stru_clustering(
        stru_esm: StructureEnsemble, 
        region_pattern: str,
        method: str = "kmeans_rmsd",
        ignore_solvent: bool = True,
        **kwargs,
    ) -> Dict[Structure, Tuple[float, float]]:
    """Perform clustering on the structure ensemble {stru_esm} and yield representative structures
    and their corresponding weigh. 

    Args:
        stru_esm: 
            A conformational ensemble of a structure.
        region_pattern: 
            A pymol-style selection pattern that defines atoms that are considered in structure-based clustering.
        method:
            Choose the algorithm of the clustering. Current supported options are:
            - kmeans_rmsd
                K-means clustering based on the RMSD value of the region selected by {region_pattern}
        ignore_solvent:
            control of solvent are removed before the calculation. True by default to speed up the calculation
            if solvent is the interest of study, set this to False.
        
        (when method = "kmeans_rmsd")
        n_clusters: 
            the desired number of clusters from K-Means.
        random_state:
            the random state of K-Means training.

    Returns:
        A dictionary that map representative structures (enzy_htp.Structure) to 
            - their weights
            - their cluster member indexes
        , respectively.
    """
    return CLUSTERING_METHOD_MAPPER[method](
        stru_esm=stru_esm,
        region_pattern=region_pattern,
        ignore_solvent=ignore_solvent,
        **kwargs,
    )

def kmeans_rmsd_clustering(
        stru_esm: StructureEnsemble, 
        region_pattern: str,
        n_clusters: int, 
        ignore_solvent: bool = True,
        random_state: int = None,
        **kwargs,
    ) -> Dict[Structure, Tuple[float, float]]:
    """K-means clustering based on the RMSD value of the region {region}"""
    result = {}
    rmsd_values = rmsd(
        stru_esm, region_pattern, ignore_solvent
    )
    X = np.array(rmsd_values).reshape(-1, 1)
    km = KMeans(n_clusters=n_clusters, random_state=random_state)
    labels = km.fit_predict(X) # label every stru in the ensemble with the cluster id

    # calculate weights of each cluster and find medoid structure
    total_num = len(labels)
    structures = list(stru_esm.structures(remove_solvent=ignore_solvent))
    for k in range(n_clusters):
        members = np.where(labels == k)[0]
        center  = km.cluster_centers_[k]
        d = np.linalg.norm(X[members] - center, axis=1)
        medoid_idx = members[d.argmin()]
        weight = len(members)/total_num

        rep_stru = structures[medoid_idx]
        result[rep_stru] = (weight, members)

    return result


CLUSTERING_METHOD_MAPPER: Dict[str, Callable[..., Dict[Structure, Tuple[float, float]]]] = {
    "kmeans_rmsd" : kmeans_rmsd_clustering
}
"""This script use EnzyHTP to calculate features from MD trajectories."""
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import multiprocessing
import logging
from pathlib import Path
import numpy as np
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.discriminant_analysis import StandardScaler
import psutil

from enzy_htp import interface as eh_interface
from enzy_htp.analysis import rmsf, coord_covariance
from enzy_htp.core.general import load_obj, save_obj
from enzy_htp.core.logger import _LOGGER as eh_logger
from enzy_htp import config as eh_config
import enzy_htp.core.file_system as fs

def feature_pca(feature_matrix: np.ndarray, n_comp: int) -> np.ndarray:
    """run PCA for a given set of features. reduce the amount of dimensions"""
    scaler = StandardScaler()
    feature_matrix_scaled = scaler.fit_transform(feature_matrix)
    pca = PCA(n_components=n_comp)
    feature_matrix_pca = pca.fit_transform(feature_matrix_scaled)

    return feature_matrix_pca

def esm_pca():
    """perform pca on esm embedding"""
    for esm_file in [
        "data/training_set/esm_features_train.pickle",
        "data/test_set/esm_features_test.pickle"
        ]:
        esm_mapper = load_obj(esm_file)
        id_list = []
        esm_matrix = []
        for i,j in esm_mapper.items():
            id_list.append(i)
            esm_matrix.append(j)
        pca_result = feature_pca(np.array(esm_matrix), 60)
        result = dict(zip(id_list, pca_result))
        print(result)
        save_obj(result, f"{Path(esm_file).with_suffix('.pca.pickle')}")

def calculate_rmsf_from_traj(traj_file: str, top_file: str, ref_pdb: str):
    """calculate RMSF for each residues from a trajectorie"""
    structure_ensemble = eh_interface.amber.load_traj(
        prmtop_path=top_file,
        traj_path=traj_file,
        ref_pdb=ref_pdb,
    )
    result = rmsf(stru_esm=structure_ensemble, by_residue=True)
    fs.safe_rm(structure_ensemble.topology_source_file)
    return list(result.values())

def unit_task_rmsf(args):
    """unit task of the parallel run"""
    idx, top_file, traj_file, ref_pdb, group, scratch = args

    eh_config.system.SCRATCH_DIR = scratch # avoid non-lock files to conflict

    result = {}
    result_file = group / "rmsf.pickle"
    print(f"working on {idx}")
    if traj_file.exists():
        if not result_file.exists():
            result[int(idx)] = calculate_rmsf_from_traj(
                top_file = str(top_file), traj_file = str(traj_file), ref_pdb = str(ref_pdb))
            print(f"finished the calculation for {idx}")
            save_obj(result, str(result_file))
        else:
            print(f"{result_file} already exists")
    else:
        print(f"MD failed in {idx}")

def create_rmsf_dataset():
    """calculate RMSF for each residues from trajectories in the dataset"""

    task_args = []
    eh_logger.setLevel(logging.WARNING)

    for master_dir in [
        "data/training_set/training_set_md/",
        "data/test_set/test_set_md/",
    ]:
        groups = Path(master_dir).glob("group_*")
        for group in groups:
            top_file = next(group.glob("traj_50/*prmtop"))
            traj_file = group / "traj_50/MD/prod.mdcrd"
            ref_pdb = next(group.glob("traj_50/*_aH.pdb"))
            idx = ref_pdb.stem.removesuffix("_rmW_rmL_rmH_aH")
            scratch = group  / "scratch"
            fs.safe_mkdir(str(scratch))
            task_args.append((idx, top_file, traj_file, ref_pdb, group, scratch))

    with multiprocessing.Pool(processes=40) as pool:
        pool.map(unit_task_rmsf, task_args)

def collect_feature_dataset(feature_file_name: str):
    """collect rmsf data from each working directory and merge to a file"""
    for master_dir in [
        "data/training_set/training_set_md/",
        "data/test_set/test_set_md/",
    ]:
        result = {}
        groups = Path(master_dir).glob("group_*")
        for group in groups:
            rmsf_file = group / feature_file_name
            if rmsf_file.exists():
                rmsf_result = load_obj(rmsf_file)
                result.update(rmsf_result)
        save_obj(result, f"{Path(master_dir).parent}/{feature_file_name}")

def rmsf_standardize():
    """How to standardize the RMSF vector?
    Method:
    1. Standardize sequence regions by dividing the sequence into #-fold
        calculate RMSF in each fold as RMS of each residue
        (RMSF_f = sqrt(mean(RMSF_i^2)))
    2. Mask-aware PCA to the same deminsion by
        Standardize the distribution of each RMSF
        zero-pad
        weighted SVD
    """
    regional_rmsf()
    svd_rmsf()

def regional_rmsf(
        rmsf_paths = (
            "data/training_set/rmsf.pickle",
            "data/test_set/rmsf.pickle",
        ),        
        n_comp: int = 20
    ):
    """
    Divide variable length RMSF by relative positions into n_comp segments
    Take RMS of each segment yield a fixed length vector
    """
    for fp in rmsf_paths:
        fp = Path(fp)
        rmsf_mapper = load_obj(str(fp)) # {id: rmsf_vector}

        out_map = {}
        for seq_idx, vec in rmsf_mapper.items():
            idx = np.linspace(0, len(vec), n_comp + 1, dtype=int) # dividing position index
            seg_means = [np.sqrt((np.array(vec[idx[i]:idx[i+1]])**2).mean()) for i in range(n_comp)]
            out_map[seq_idx] = np.array(seg_means, dtype=float)

        # save
        out_fp = fp.with_stem(fp.stem + f"_seg{n_comp}")
        save_obj(out_map, out_fp)

def svd_rmsf(
        rmsf_paths = (
            "data/training_set/rmsf.pickle",
            "data/test_set/rmsf.pickle",
        ),        
        n_comp: int = 60
    ):
    """perform weighted SVD on esm embedding.
    Set the weight of masks to be 0
    NOTE this will lose the overall RMSD difference btw always concat with other features"""
    for rmsf_file in rmsf_paths:
        rmsf_mapper = load_obj(rmsf_file)
        id_list, rmsf_matrix = zip(*rmsf_mapper.items())
        l_max = max(map(len, rmsf_matrix))
        n_seq = len(rmsf_matrix)

        # zero-pad
        X_pad  = np.zeros((n_seq, l_max), dtype=float)
        Mask   = np.zeros_like(X_pad, bool)        
        for i, vec in enumerate(rmsf_matrix):
            vec = np.array(vec)
            ln = len(vec)
            X_pad[i, :ln] = (vec - vec.mean()) / (vec.std(ddof=0) + 1e-8)  # z-score
            Mask [i, :ln] = True

        # Mask-aware centralize
        col_sum   = (X_pad * Mask).sum(axis=0, keepdims=True)
        col_count = Mask.sum(axis=0, keepdims=True).clip(min=1)
        X_center  = X_pad - col_sum / col_count

        # SVD rmsf_matrix to 60 components
        X_weighted = X_center * np.sqrt(Mask)
        svd = TruncatedSVD(n_components=n_comp, algorithm="randomized", random_state=42)
        Z   = svd.fit_transform(X_weighted)
        result = {i: z for i, z in zip(id_list, Z)}
        save_obj(result, f"{Path(rmsf_file).with_suffix('.svd60.pickle')}")

def create_covariance_dataset():
    """calculate covariance matrix from trajectories in the dataset"""

    task_args = []
    eh_logger.setLevel(logging.WARNING)

    for master_dir in [
        "data/training_set/training_set_md/",
        "data/test_set/test_set_md/",
    ]:
        groups = Path(master_dir).glob("group_*")
        for group in groups:
            top_file = next(group.glob("traj_50/*prmtop"))
            traj_file = group / "traj_50/MD/prod.mdcrd"
            ref_pdb = next(group.glob("traj_50/*_aH.pdb"))
            idx = ref_pdb.stem.removesuffix("_rmW_rmL_rmH_aH")
            scratch = group  / "scratch"
            fs.safe_mkdir(str(scratch))
            task_args.append((idx, top_file, traj_file, ref_pdb, group, scratch))

    with multiprocessing.Pool(processes=40) as pool:
        pool.map(unit_task_covariance, task_args)

def calculate_covariance_from_traj(traj_file: str, top_file: str, ref_pdb: str) -> np.ndarray:
    """calculate covariance matrix from a trajectorie"""
    structure_ensemble = eh_interface.amber.load_traj(
        prmtop_path=top_file,
        traj_path=traj_file,
        ref_pdb=ref_pdb,
    )
    result = coord_covariance(
        stru_esm=structure_ensemble,
        region_pattern="polymer and (n. CA)",
        reference_type = "average",
        mass_weighted = False,
    )
    fs.safe_rm(structure_ensemble.topology_source_file)
    return result

def unit_task_covariance(args):
    """unit task of the parallel run"""
    idx, top_file, traj_file, ref_pdb, group, scratch = args
    proc = psutil.Process()

    eh_config.system.SCRATCH_DIR = scratch # avoid non-lock files to conflict

    result = {}
    result_file = group / "covariance.pickle"
    print(f"working on {idx}")
    if traj_file.exists():
        if not result_file.exists():
            result[int(idx)] = calculate_covariance_from_traj(
                top_file = str(top_file), traj_file = str(traj_file), ref_pdb = str(ref_pdb))
            print(f"finished the calculation for {idx}")
            save_obj(result, str(result_file))
        else:
            print(f"{result_file} already exists")
    else:
        print(f"MD failed in {idx}")
    print(f"[{idx}] num_thread used: {proc.num_threads()}")

def covariance_standardize():
    """How to standardize the covariance matrix?
    Method:
    1. Calculate the eigen-spectrum and keep the largest 40 eigenvalues
    2. Create the DCCM-Histogram of the covariance matrix (24 bin)
    3. Flat the upper triangle into a vector and do weighted SVD on vectors
        in the dataset similar to RMSF. (n_comp=60)
    """
    eigen_spec_covariance()
    hist_covariance()
    svd_flat_covariance(remove_diag=False)
    svd_flat_covariance(remove_diag=True)


def _top_k_eigs(C: np.ndarray,
                k_keep: int,
                drop_rigid: int,
                log_scale: bool) -> np.ndarray:

    lamb = np.linalg.eigvalsh(np.asarray(C, float))
    lamb = lamb[drop_rigid:]
    vec  = lamb[-k_keep:][::-1]
    if len(vec) < k_keep:
        vec = np.pad(vec, (0, k_keep - len(vec)))
    if log_scale:
        vec = np.log1p(np.maximum(vec, 0.0))
    return vec

def eigen_spec_covariance(
        train_data = "data/training_set/covariance.pickle",
        test_data = "data/test_set/covariance.pickle",
        k_keep: int   = 40,
        drop_rigid: int = 6,
        log_scale: bool = True,
        suffix: str   = ".eig40.pickle",
    ):
    """
    For each covariance matrix C (3N*3N), compute eigenvalues λ_i,
    discard the smallest `drop_rigid` (≈ rigid-body modes after alignment),
    keep the next `k_keep` largest, pad with 0 if needed.
    Optionally apply log-scale and global z-score.

    Save result in a file:
        {id: 1D ndarray}
    """
    train_mapper = load_obj(train_data)           # {id: C}
    train_vecs   = {
        pid: _top_k_eigs(C, k_keep, drop_rigid, log_scale)
        for pid, C in train_mapper.items()
    }
    scaler = StandardScaler()
    train_mat_z = scaler.fit_transform(np.vstack(list(train_vecs.values())))
    for pid, z in zip(train_vecs, train_mat_z):
        train_vecs[pid] = z.astype(np.float32)

    save_obj(train_vecs, Path(train_data).with_suffix(suffix))

    # test set transform using the same set up as training set standardized
    test_mapper = load_obj(test_data)
    test_vecs = {}
    for pid, C in test_mapper.items():
        vec = _top_k_eigs(C, k_keep, drop_rigid, log_scale)
        test_vecs[pid] = scaler.transform(vec[None, :]).ravel().astype(np.float32)

    save_obj(test_vecs, Path(test_data).with_suffix(suffix))

def hist_covariance(
        data_paths = (
            "data/training_set/covariance.pickle",
            "data/test_set/covariance.pickle",
        ),
        bins       = np.linspace(-1.0, 1.0, 25),
        density    = True,
        suffix     = ".dccm24.pickle",
    ):
    """
    Convert each covariance matrix to a fixed-length DCCM histogram feature.

    For every (NxN) covariance matrix C:
      1) Normalize: R_ij = C_ij / sqrt(C_ii * C_jj)
      2) Flatten: R_ij (i < j) as 1-D vector
      3) Calculate the histogram of these vectors in bins (density=True → prob density)
    The resulting 1D array has length = len(bins) - 1  (defaul: 24).

    Save results to a file:
        {protein_id: histogram_vector}
    """
    for data_file in data_paths:
        data_mapper = load_obj(data_file)          # {id: ndarray (3N x 3N)}
        id_list, cov_mats = zip(*data_mapper.items())

        feats = {}
        for pid, C in zip(id_list, cov_mats):
            C = np.asarray(C, dtype=float)
            if C.shape[0] != C.shape[1]:
                raise ValueError(f"{pid}: covariance must be square, got {C.shape}")

            # 1. Normalize
            diag = np.diag(C).copy()
            if np.any(diag <= 0):
                diag = np.where(diag <= 0, 1e-8, diag)
            inv_sqrt = 1.0 / np.sqrt(diag)
            R = C * inv_sqrt[:, None] * inv_sqrt[None, :]

            # 2. Flat
            triu_vals = R[np.triu_indices_from(R, k=1)]

            # 3. Histogram
            hist, _ = np.histogram(triu_vals, bins=bins, density=density)
            feats[pid] = hist.astype(np.float32)

        out_file =str(Path(data_file).with_suffix(suffix))
        save_obj(feats, out_file)

def svd_flat_covariance(        
        data_paths = (
            "data/training_set/covariance.pickle",
            "data/test_set/covariance.pickle",
        ),        
        n_comp: int = 60,
        remove_diag: bool = True,
    ):
    """perform weighted SVD on flattened upper triangle of the covariance matrix.
    Set the weight of masks of shorter vectors to be 0"""
    for data_file in data_paths:
        data_mapper = load_obj(data_file)
        id_list, raw_data_matrix = zip(*data_mapper.items())
        # TODO here flat every matrix in data_matrix to a vector of upper triangle
        data_matrix = []
        for M in raw_data_matrix:
            M = np.asarray(M, dtype=np.float32)
            iu = np.triu_indices(M.shape[0], k=int(remove_diag))
            data_matrix.append(M[iu])

        l_max = max(map(len, data_matrix))
        n_seq = len(data_matrix)

        # zero-pad
        X_pad  = np.zeros((n_seq, l_max), dtype=float)
        Mask   = np.zeros_like(X_pad, bool)        
        for i, vec in enumerate(data_matrix):
            vec = np.array(vec)
            ln = len(vec)
            X_pad[i, :ln] = (vec - vec.mean()) / (vec.std(ddof=0) + 1e-8)  # z-score
            Mask [i, :ln] = True

        # Mask-aware centralize
        col_sum   = (X_pad * Mask).sum(axis=0, keepdims=True)
        col_count = Mask.sum(axis=0, keepdims=True).clip(min=1)
        X_center  = X_pad - col_sum / col_count

        # SVD data_matrix to 60 components
        X_weighted = X_center * np.sqrt(Mask)
        svd = TruncatedSVD(n_components=n_comp, algorithm="randomized", random_state=42)
        Z   = svd.fit_transform(X_weighted)
        result = {i: z for i, z in zip(id_list, Z)}

        out_file = str(Path(data_file).with_suffix(
            f'.svd{n_comp}k{int(remove_diag)}.pickle'
        ))
        
        save_obj(result, out_file)

def main():
    # esm_pca()
    # collect_rmsf_dataset()
    covariance_standardize()
    # create_covariance_dataset()
    # collect_feature_dataset("covariance.pickle")
    pass

if __name__ == "__main__":
    main()

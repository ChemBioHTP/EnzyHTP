"""Integration test for umbrella sampling workflow.

This test demonstrates a complete umbrella sampling workflow from start to finish:
1. Run umbrella sampling simulations
2. Extract reaction coordinates
3. Calculate probability densities
4. Run WHAM analysis
5. Generate plots

Author: GitHub Copilot (with guidance from EnzyHTP conventions)
Date: 2026-03-10
"""

import os
import glob
import pytest
import numpy as np
from functools import partial
from pathlib import Path

from enzy_htp.core.clusters.accre import Accre
import enzy_htp.core.file_system as fs
from enzy_htp.structure import structure_constraint as stru_cons
from enzy_htp.geometry import umbrella_sampling
from enzy_htp.analysis import (
    extract_reaction_coordinate,
    probability_density,
    wham_pmf,
    plot_probability_density,
    plot_pmf,
)
from enzy_htp import interface
from enzy_htp import PDBParser

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"
sp = PDBParser()


def _get_missing_wham_executables() -> list:
    """Get missing WHAM executables from the registered EnzyHTP interface."""
    if not hasattr(interface, "wham"):
        return ["wham interface not registered"]

    return interface.wham.missing_executables()


@pytest.mark.accre
@pytest.mark.long
@pytest.mark.integration
def test_umbrella_sampling_full_workflow():
    """Complete umbrella sampling workflow integration test.
    
    This test runs:
    1. Umbrella sampling with 3 windows
    2. Extracts reaction coordinates from trajectories
    3. Calculates probability densities
    4. Runs WHAM to get PMF (if WHAM is available)
    5. Generates plots
    
    This is a comprehensive test that may take several minutes to run.
    """
    # Setup
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J": (0, 1)})
    
    test_param_method = interface.amber.build_md_parameterizer(
        ncaa_param_lib_path=f"{WORK_DIR}/ncaa_lib",
    )
    
    cluster_job_config = {
        "cluster": Accre(),
        "res_keywords": {
            "account": "csb_gpu_acc",
            "partition": "turing",
        }
    }
    
    # Force constant for umbrella potential (kcal/mol/Å²)
    force_constant = 5.0
    
    # Define constraint generator
    constraint_gen = partial(
        stru_cons.create_group_distance_constraint,
        "resi 254",  # Ligand selection
        "resi 9+11+48+50+101+128+169+201+202+222+224 & n. C+CA+N",  # Enzyme selection
        params={
            "amber": {
                "ialtd": 0,
                "r1": "x-999999",
                "r2": "x",
                "r3": "x",
                "r4": "x+999999",
                "rk2": force_constant,
                "rk3": force_constant,
                "rs_filepath": "{mdstep_dir}/0.rs",
            }
        }
    )

    cv = AmberConstraintCollectiveVariable(constraint_gen=constraint_gen)

    
    # ===== Step 1: Run Umbrella Sampling =====
    # window_targets = [4.0, 5.0, 6.0]  # Three windows (Å)
    window_targets = cv.generate_window_targets(
        [4.0, 5.0, 6.0]
    )
    window_targets = cv.generate_window_targets(
        start=4.0,
        end=6.0,
        step=1.0,
    )

    umbrella_work_dir = f"{WORK_DIR}integration_umbrella/"
    
    umbrella_result: UmbrellaSamplingResult = umbrella_sampling(
        stru=test_stru,
        param_method=test_param_method,
        cv=cv,
        window_targets=window_targets,
        parallel_runs=1,
        work_dir=umbrella_work_dir,
        prod_time=0.05,  # Short for testing (50 ps)
        record_period=0.005,  # Record every 5 ps
        cluster_job_config=cluster_job_config,
        job_check_period=10,
    )
    
    assert len(umbrella_result.windows) == 3, "Should have 3 windows"
    assert len(umbrella_result.final_structures) == 3, "Should have 3 final structures"
    
    # ===== Step 2: Extract Reaction Coordinates =====
    rc_by_window = []
    timeseries_files = []
    
    for window_idx in range(len(window_targets)):
        # Get the  first replica ensemble for this window
        ensemble = umbrella_result.windows[window_idx].replicas[0]
        
        # Extract RC timeseries
        rc_values = extract_reaction_coordinate(
            ensemble=ensemble,
            cv=cv,
            engine="cpptraj"
        )
        
        rc_by_window.append(rc_values)
        
        # Save timeseries to file for WHAM
        # ==> modularize
        ts_file = f"{umbrella_work_dir}window_{window_idx}_timeseries.dat"
        with open(ts_file, 'w') as f:
            f.write("#Frame Distance\n")
            for frame_idx, value in enumerate(rc_values):
                f.write(f"{frame_idx+1} {value:.6f}\n")
        timeseries_files.append(ts_file)
        
        print(f"Window {window_idx}: extracted {len(rc_values)} frames, "
              f"mean RC = {rc_values.mean():.2f} Å")
    
    # ===== Step 3: Calculate Probability Densities =====
    # ==> modularize
    all_densities = []
    for window_idx, rc_values in enumerate(rc_by_window):
        centers, density = probability_density(rc_values, bins=30)
        all_densities.append((centers, density))
        print(f"Window {window_idx}: peak at {centers[np.argmax(density)]:.2f} Å")
    
    # ===== Step 5: Run WHAM (if available) =====
    missing_wham_exes = _get_missing_wham_executables()

    if not missing_wham_exes:
        print("WHAM executable found, running PMF calculation...")
        
        rc, pmf = wham_pmf(
            timeseries_files=timeseries_files,
            target_distances=window_targets,
            force_constant=force_constant,
            temperature=300.0,
            bins=50,
            tolerance=0.0001,
            work_dir=f"{umbrella_work_dir}wham_analysis/",
        )
        
        print(f"PMF calculated: min = {pmf.min():.2f}, max = {pmf.max():.2f} kcal/mol")
        
        
        # Verify PMF properties
        assert pmf.min() == pytest.approx(0.0, abs=0.01), "PMF minimum should be ~0"
        assert pmf.max() >= 0.0, "PMF should be non-negative"
        assert 8 <= pmf.max() <= 11.0, "PMF maximum should be in expected range (8-11 kcal/mol)"
        assert len(rc) == 50, "Should have 50 bins"
        
    else:
        print("WHAM executable not found in registered interface, skipping PMF calculation")
        pytest.skip(f"WHAM not available for PMF calculation: {missing_wham_exes}")
    
    # ===== Cleanup =====
    fs.safe_rmdir(umbrella_work_dir)
    fs.clean_temp_file_n_dir(
        glob.glob("slurm-*.out") +
        glob.glob(f"{WORK_DIR}/ncaa_lib/H5J*")
    )


@pytest.mark.integration
def test_umbrella_sampling_analysis_only():
    """Test analysis functions with pre-generated synthetic data.
    
    This test doesn't require running actual MD simulations or having
    WHAM installed. It uses synthetic data to test the analysis pipeline.
    """
    work_dir = Path(WORK_DIR) / "synthetic_umbrella"
    work_dir.mkdir(parents=True, exist_ok=True)
    
    # ===== Generate synthetic umbrella sampling data =====
    # Simulate 3 windows with gaussian distributions
    window_targets = [4.0, 5.5, 7.0]
    force_constant = 10.0
    
    rc_by_window = []
    timeseries_files = []
    
    for window_idx, target in enumerate(window_targets):
        # Generate synthetic data: gaussian centered at target with some width
        rc_values = np.random.normal(loc=target, scale=0.3, size=500)
        rc_by_window.append(rc_values)
        
        # Save to file
        ts_file = work_dir / f"window_{window_idx}.dat"
        with open(ts_file, 'w') as f:
            f.write("#Frame Distance\n")
            for frame_idx, value in enumerate(rc_values):
                f.write(f"{frame_idx+1} {value:.6f}\n")
        timeseries_files.append(str(ts_file))
    
    # ===== Test probability density calculation =====
    centers, combined_density = probability_density(rc_by_window, bins=50)
    assert len(centers) == 50
    assert len(combined_density) == 50
    
    # Test individual densities
    for rc_values in rc_by_window:
        centers, density = probability_density(rc_values, bins=30)
        assert len(centers) == 30
        # Integral should be approximately 1
        bin_width = centers[1] - centers[0]
        integral = np.sum(density) * bin_width
        assert np.isclose(integral, 1.0, atol=0.1)
    
    # ===== Test plotting =====
    prob_plot = work_dir / "prob_density_synthetic.png"
    plot_probability_density(
        rc_values_by_window=rc_by_window,
        output_path=str(prob_plot),
        bins=30,
    )
    assert prob_plot.exists()
    
    # ===== Test WHAM if available =====
    missing_wham_exes = _get_missing_wham_executables()

    if not missing_wham_exes:
        rc, pmf = wham_pmf(
            timeseries_files=timeseries_files,
            target_distances=window_targets,
            force_constant=force_constant,
            temperature=300.0,
            bins=40,
            work_dir=str(work_dir / "wham"),
        )
        
        # Plot PMF
        pmf_plot = work_dir / "pmf_synthetic.png"
        plot_pmf(
            reaction_coordinate=rc,
            pmf=pmf,
            output_path=str(pmf_plot),
        )
        assert pmf_plot.exists()
        
        # Verify PMF is reasonable
        assert len(rc) == 40
        assert len(pmf) == 40
        assert pmf.min() == pytest.approx(0.0, abs=0.01)
    else:
        pytest.skip(f"WHAM not available for PMF calculation: {missing_wham_exes}")
    
    # ===== Cleanup =====
    fs.safe_rmdir(str(work_dir))

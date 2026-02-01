"""Usage: python af2_main.py <cpu_model_name>

Utilities:
- ``collect_sequence_from_uniprot``: read ``apo_pdb_uniprot.txt``, fetch AlphaFold
  prediction metadata from the public API, and save a CSV plus ``input.fasta`` in
  the same order as the input list.
- ``run_af2_main``: submit AlphaFold2 jobs using enzy_htp.
"""
import csv
import json
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Optional, Tuple
import pandas as pd

from enzy_htp.structure_prediction import predict_structure
from enzy_htp.core.clusters.accre_r9 import AccreR9
from enzy_htp.core.job_manager import ClusterJobConfig
from enzy_htp import config as eh_config


ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
INPUT_TABLE = Path("apo_pdb_uniprot.txt")
OUTPUT_CSV = Path("uniprot_sequences.csv")
OUTPUT_FASTA = Path("input.fasta")


def _seq_missing_files(seq_dir: Path) -> list[str]:
    """Return the missing unrelaxed PDBs for a sequence directory."""
    required = [f"unrelaxed_model_{i}_ptm_pred_0.pdb" for i in range(1, 6)]
    return [fname for fname in required if not (seq_dir / fname).exists()]


def _fetch_alphafold_record(uniprot_id: str) -> Optional[dict]:
    """Return the best AlphaFold record for a UniProt ID or ``None`` on failure."""

    url = ALPHAFOLD_API.format(uniprot_id=uniprot_id)
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            payload = json.loads(resp.read())
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        print(f"[WARN] Failed to fetch {uniprot_id}: {exc}", file=sys.stderr)
        return None

    if not isinstance(payload, list) or not payload:
        print(f"[WARN] Unexpected response for {uniprot_id}: {payload}", file=sys.stderr)
        return None

    # Prefer the entry with the highest version; if tied, prefer the F1 isoform.
    def _score(entry: dict) -> Tuple[int, int]:
        version = entry.get("latestVersion", 0) or 0
        is_f1 = 1 if str(entry.get("modelEntityId", "")).endswith("F1") else 0
        return version, is_f1

    return max(payload, key=_score)


def collect_sequence_from_uniprot(
    input_table: Path = INPUT_TABLE,
    output_csv: Path = OUTPUT_CSV,
    output_fasta: Path = OUTPUT_FASTA,
    ref_table: Optional[Path] = None,
):
    """Fetch sequences/tool metadata from AlphaFold and save CSV + FASTA.

    The CSV columns are: ``pdb_id``, ``uniprot_id``, ``toolUsed``, ``msaUrl``,
    ``sequence``. Rows keep the same order as ``apo_pdb_uniprot.txt``.
    The FASTA file writes only entries with a retrieved sequence.
    """

    rows = []
    seen_pairs = set()
    fasta_lines: list[str] = []

    if ref_table:
        ref_table = pd.read_csv(ref_table)
        seen_pairs = set(zip(ref_table["pdb_id"], ref_table["uniprot_id"]))

    with input_table.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if stripped.lower().startswith("pdb_id"):
                # header line
                continue

            parts = stripped.split()
            if len(parts) < 2:
                print(f"[WARN] Skip malformed line: {line.rstrip()}", file=sys.stderr)
                continue

            pdb_id, uniprot_id = parts[0], parts[1]
            key = (pdb_id, uniprot_id)
            if key in seen_pairs:
                print(f"[INFO] Skip duplicate entry: {pdb_id} {uniprot_id}", file=sys.stderr)
                continue
            seen_pairs.add(key)

            record = _fetch_alphafold_record(uniprot_id)

            sequence = record.get("sequence") if record else ""
            tool_used = record.get("toolUsed") if record else ""
            msa_url = record.get("msaUrl") if record else ""

            rows.append(
                {
                    "pdb_id": pdb_id,
                    "uniprot_id": uniprot_id,
                    "toolUsed": tool_used,
                    "msaUrl": msa_url,
                    "sequence": sequence,
                }
            )

            if sequence:
                fasta_lines.append(f">{pdb_id}|{uniprot_id}")
                fasta_lines.append(sequence)

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="") as csvfile:
        writer = csv.DictWriter(
            csvfile, fieldnames=["pdb_id", "uniprot_id", "toolUsed", "msaUrl", "sequence"]
        )
        writer.writeheader()
        writer.writerows(rows)

    if fasta_lines:
        with output_fasta.open("w") as fasta_file:
            fasta_file.write("\n".join(fasta_lines) + "\n")

    print(f"[INFO] Wrote {len(rows)} rows to {output_csv} and FASTA to {output_fasta}")


def count_completed_sequences(work_dir: Path) -> tuple[int, int]:
    """Count how many sequence jobs have produced all 5 unrelaxed model PDBs.

    A job is considered complete if the directory ``seq_*`` inside ``work_dir``
    contains the five files ``unrelaxed_model_<1-5>_ptm_pred_0.pdb``.
    Returns a tuple of ``(completed, total)`` and prints a short summary plus
    the missing files (if any).
    """

    work_dir = Path(work_dir)
    seq_dirs = [p for p in work_dir.glob("seq_*") if p.is_dir()]

    completed = 0
    incomplete: list[tuple[str, list[str]]] = []
    for seq_dir in seq_dirs:
        missing = _seq_missing_files(seq_dir)
        if missing:
            incomplete.append((seq_dir.name, missing))
        else:
            completed += 1

    print(f"[INFO] Completed sequences: {completed}/{len(seq_dirs)} in {work_dir}")
    if incomplete:
        print("[INFO] Incomplete sequence directories (missing files):")
        for name, missing in incomplete:
            print(f"  {name}: {', '.join(missing)}")
    return completed, len(seq_dirs)


def _parse_sequences_from_submit(script_path: Path) -> list[int]:
    """Extract sequence indices from a submit_alphafold_X.cmd file."""

    import re

    try:
        text = script_path.read_text(errors="ignore")
    except FileNotFoundError:
        return []

    return sorted(set(int(m) for m in re.findall(r"seq_(\d+)\.fasta", text)))


def _error_snippet_from_log(log_text: str) -> Optional[str]:
    """Return a short error snippet if a Slurm log indicates failure."""

    markers = [
        "Traceback (most recent call last)",
        "RuntimeError",
        "CUDA_ERROR",
        "Killed",
        "ERROR",
        "Error",
        "exception",
        "out of memory",
    ]
    lower = log_text.lower()
    needle = None
    for m in markers:
        if m.lower() in lower:
            needle = m
            break

    if not needle:
        return None

    lines = log_text.splitlines()
    idx = max(i for i, ln in enumerate(lines) if needle.lower() in ln.lower())
    snippet = lines[idx : min(len(lines), idx + 8)]
    return "\n".join(snippet)

def analyze_submitted_jobs(
    submitted_log: Path = Path("submitted_job_ids.log"),
    log_dir: Path = Path("."),
) -> None:
    """Analyze Slurm jobs listed in submitted_job_ids.log.

    For each job id and submit script:
    - find the corresponding slurm-<jobid>.out
    - infer seq indices from the submit script
    - check required outputs (5 unrelaxed PDBs) under the script's parent dir
    - print status and reason for failures
    """

    if not submitted_log.exists():
        print(f"[ERROR] Missing submitted job list: {submitted_log}", file=sys.stderr)
        sys.exit(1)

    rows = []
    for line in submitted_log.read_text().splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split()
        if len(parts) < 2:
            continue
        job_id, script_rel = parts[0], parts[1]
        rows.append((job_id, Path(script_rel)))

    if not rows:
        print(f"[WARN] No jobs found in {submitted_log}")
        return

    ok = fail = missing_log = running = pending = 0
    for job_id, script_rel in rows:
        script_path = script_rel if script_rel.is_absolute() else submitted_log.parent / script_rel
        work_dir = script_path.parent
        seq_ids = _parse_sequences_from_submit(script_path)

        log_path = log_dir / f"slurm-{job_id}.out"

        log_text = log_path.read_text(errors="ignore")
        err_snippet = _error_snippet_from_log(log_text) or ""

        missing_outputs: list[int] = []
        for seq_id in seq_ids:
            seq_dir = work_dir / f"seq_{seq_id}"
            if not seq_dir.exists() or _seq_missing_files(seq_dir):
                missing_outputs.append(seq_id)

        state, state_kw = AccreR9.get_job_state(job_id)

        if not err_snippet and not missing_outputs:
            ok += 1
            if state and state != "complete":
                # Outputs look fine but job still running/pending? note it.
                print(f"[ OK?] {job_id} ({script_rel}) | seqs: {seq_ids or '-'} | state={state}")
            else:
                print(f"[ OK ] {job_id} ({script_rel}) | seqs: {seq_ids or '-'}")
        elif state in {"run", "pend"}:
            if state == "run":
                running += 1
            else:
                pending += 1
            print(f"[{state.upper():5}] {job_id} ({script_rel}) | seqs: {seq_ids or '-'} | "
                  f"missing_outputs: {missing_outputs or '-'}")
        else:
            fail += 1
            print(f"[FAIL] {job_id} ({script_rel}) | seqs: {seq_ids or '-'} | "
                  f"missing_outputs: {missing_outputs or '-'} | state={state if state != 'unknown' else state_kw}")
            if err_snippet:
                print("       Error snippet:")
                for ln in err_snippet.splitlines():
                    print(f"       {ln}")

    total = ok + fail + missing_log + running + pending
    print(f"\n[SUMMARY] ok={ok}, fail={fail}, running={running}, pending={pending}, no_log={missing_log}, total_jobs={total}")


def generate_rerun_scripts(
    submitted_log: Path = Path("submitted_job_ids.log"),
    out_dir_name: str = "rerun",
    mem: str = "60G",
) -> None:
    """Generate rerun submit scripts with completed seqs removed and higher memory.

    For each script listed in ``submitted_log``:
    - detect incomplete ``seq_<id>`` under the script's work directory
    - rewrite the python line ``--fasta_paths`` to include only incomplete seqs
    - bump ``#SBATCH --mem`` to ``mem``
    - write a sibling script named ``<orig>_rerun.cmd``
    Scripts with all seqs complete are skipped.
    """

    if not submitted_log.exists():
        print(f"[ERROR] Missing submitted job list: {submitted_log}", file=sys.stderr)
        sys.exit(1)

    rows: list[Path] = []
    for line in submitted_log.read_text().splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split()
        if len(parts) < 2:
            continue
        script_rel = Path(parts[1])
        job_id = parts[0]
        rows.append((job_id, script_rel) if script_rel.is_absolute() else (job_id, submitted_log.parent / script_rel))

    if not rows:
        print(f"[WARN] No scripts found in {submitted_log}")
        return

    def _python_line_indices(lines: list[str]) -> list[int]:
        return [i for i, ln in enumerate(lines) if "--fasta_paths" in ln and ln.lstrip().startswith("python")]

    for job_id, script_path in rows:
        if not script_path.exists():
            print(f"[WARN] Script not found, skip: {script_path}")
            continue

        seq_ids = _parse_sequences_from_submit(script_path)
        if not seq_ids:
            print(f"[WARN] No seq ids detected in {script_path}, skip")
            continue

        work_dir = script_path.parent
        out_dir = work_dir / out_dir_name
        out_dir.mkdir(exist_ok=True)
        incomplete_ids: list[int] = []
        for seq_id in seq_ids:
            seq_dir = work_dir / f"seq_{seq_id}"
            if seq_dir.exists() and not _seq_missing_files(seq_dir):
                continue  # already complete
            incomplete_ids.append(seq_id)

        if not incomplete_ids:
            print(f"[SKIP] all seqs complete for {script_path}")
            continue

        job_state, job_state_kw = AccreR9.get_job_state(job_id)
        if job_state in {"run", "pend"}:
            print(f"[SKIP] job {job_id} still {job_state}, skip rerun script generation")
            continue

        lines = script_path.read_text().splitlines()

        # update memory line
        new_lines: list[str] = []
        for ln in lines:
            if ln.startswith("#SBATCH --mem="):
                new_lines.append(f"#SBATCH --mem={mem}")
            else:
                new_lines.append(ln)

        py_idxs = _python_line_indices(new_lines)
        if not py_idxs:
            print(f"[WARN] No python line with --fasta_paths in {script_path}, skip")
            continue

        # assume first match is the target command
        idx = py_idxs[0]
        py_line = new_lines[idx]

        parts = py_line.split("--fasta_paths", 1)
        if len(parts) < 2:
            print(f"[WARN] malformed python line in {script_path}, skip")
            continue

        left = parts[0] + "--fasta_paths "
        rhs = parts[1]
        rhs_parts = rhs.split(None, 1)  # value, rest
        tail = rhs_parts[1] if len(rhs_parts) > 1 else ""

        new_fasta = ",".join(f"{work_dir.name}/seq_{sid}.fasta" for sid in incomplete_ids)
        rebuilt = left + new_fasta + (" " + tail if tail else "")
        new_lines[idx] = rebuilt

        out_path = out_dir / (script_path.stem + "_rerun.cmd")
        out_path.write_text("\n".join(new_lines) + "\n")
        print(f"[WRITE] {out_path} | kept seqs {incomplete_ids}")

def run_af2_main():
    """run AF2 structure prediction on Accre.
    Usage: python af2_main.py <cpu_model_name> <uniprot_sequences.csv>
    Example: python af2_main.py nvidia_rtx_a6000 part1.csv
    """
    seq_table = pd.read_csv(sys.argv[2])  # e.g. uniprot_sequences.csv
    seq_list = seq_table["sequence"].tolist()

    cluster_job_config = ClusterJobConfig(
        cluster = AccreR9(),
        res_keywords = {
            "account": "csb_gpu_acc",
            "partition": "batch_gpu",
            "node_cores": f"{sys.argv[1]}:1",
            "walltime": "3-00:00:00",
            "mem_per_core": "60G",
            # "account": "yang_lab_csb_iacc",
            # "partition": "interactive_gpu",
            # "qos": "debug_iacc",
            # "node_cores": "nvidia_rtx_a4000:1",
            # "walltime": "30:00",
            }
        )

    eh_config.alphafold.INSTALL_TYPE="alphafold2_native_python"
    eh_config.alphafold.EXECUTABLE_PATH="/sb/apps/alphafold232/alphafold/run_alphafold.py"
    eh_config.alphafold.DATA_DIR="/sb/apps/alphafold-data.230"
    work_dir = sys.argv[1] + "_af2_results"

    predict_structure(
        sequences=seq_list, 
        cluster_job_config=cluster_job_config,
        seq_per_job=6,
        array_size=100,
        num_relax=0,
        # use_precomputed_msas=True, # TODO support automatically copy files to work_dir/msas/ 
        # precomputed_msa_path_list = ["xxx.a3m"]
        work_dir=work_dir,
        use_templates=True,
        max_template_date="2021-02-15",
        )

def main():
    # collect_sequence_from_uniprot(input_table=Path("obsolete_uniprot_sequences.txt"),
    #                              output_csv=Path("obsolete_uniprot_sequences.csv"),
    #                              ref_table="filtered_uniprot_sequences.csv")
    if len(sys.argv) >= 2 and sys.argv[1] in {"count", "--count", "count_completed"}:
        if len(sys.argv) < 3:
            print("Usage: python af2_main.py count <work_dir>", file=sys.stderr)
            sys.exit(1)
        count_completed_sequences(Path(sys.argv[2]))
        return

    if len(sys.argv) >= 2 and sys.argv[1] in {"analyze", "analyze_jobs"}:
        submitted = Path(sys.argv[2]) if len(sys.argv) >= 3 else Path("submitted_job_ids.log")
        log_dir = Path(sys.argv[3]) if len(sys.argv) >= 4 else Path(".")
        analyze_submitted_jobs(submitted, log_dir)
        return

    # if len(sys.argv) < 3:
    #     print("Usage: python af2_main.py <cpu_model_name> <uniprot_sequences.csv>", file=sys.stderr)
    #     print("Or:    python af2_main.py count <work_dir>")
    #     print("Or:    python af2_main.py analyze [submitted_job_ids.log] [log_dir]")
    #     sys.exit(1)
    generate_rerun_scripts()
    # run_af2_main()

if __name__ == "__main__":
    main()

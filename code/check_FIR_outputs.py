#!/usr/bin/env python3
"""
QC: Verify that FIR L1 outputs exist (stats/cope1.nii.gz) for BOTH tasks (mid, sharedreward),
while honoring an explicit (sub, ses) skip list (mirrors your qsub driver).

Assumptions:
- Run from the code directory (script lives there).
- Project root is the parent of the code directory.
- Expected runs are inferred from derivatives/fsl/EVFiles-FIR/...
- A cope1 is only "expected" when:
    (a) EV run dir exists, and
    (b) the corresponding fMRIPrep preproc BOLD exists, and
    (c) the corresponding confounds TSV exists, and
    (d) the session is NOT in SKIP.
Outputs: prints only summaries + what’s missing.
"""

from __future__ import annotations

import re
from pathlib import Path
from collections import defaultdict


# --------------------------- Config ---------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR   = SCRIPT_DIR.parent
FSL_DERIV  = ROOT_DIR / "derivatives" / "fsl"
EV_BASE    = FSL_DERIV / "EVFiles-FIR"
FMRIPREP   = ROOT_DIR / "derivatives" / "fmriprep"
CNF_DIR    = FSL_DERIV / "confounds_tedana"

TASKS      = ["mid", "sharedreward"]
ECHOES     = ["single-echo", "multi-echo"]
CNFDS      = ["cnfds-fmriprep", "cnfds-tedana"]

# Explicit skip list (sub, ses) pairs — matches your bash snippet
# Format: ("101","04") corresponds to sub-101 ses-04
SKIP = {
    ("101", "04"),
    ("101", "05"),
    ("101", "12"),
    ("103", "12"),
}


def _parse_id(name: str, prefix: str) -> str | None:
    m = re.match(rf"^{re.escape(prefix)}-(.+)$", name)
    return m.group(1) if m else None


def _bold_path(sub: str, ses: str, task: str, run: str, echo: str) -> Path:
    func = FMRIPREP / f"sub-{sub}" / f"ses-{ses}" / "func"
    if echo == "multi-echo":
        return func / f"sub-{sub}_ses-{ses}_task-{task}_run-{run}_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz"
    if echo == "single-echo":
        return func / f"sub-{sub}_ses-{ses}_task-{task}_run-{run}_echo-2_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz"
    raise ValueError(f"Unexpected echo: {echo}")


def _confounds_path(sub: str, ses: str, task: str, run: str, cnfd: str) -> Path:
    base = CNF_DIR / f"sub-{sub}" / f"ses-{ses}"
    if cnfd == "cnfds-tedana":
        return base / f"sub-{sub}_ses-{ses}_task-{task}_run-{run}_desc-TedanaPlusConfounds.tsv"
    return base / f"sub-{sub}_ses-{ses}_task-{task}_run-{run}_desc-fslConfounds.tsv"


def _feat_dir(sub: str, ses: str, task: str, run: str, echo: str, cnfd: str) -> Path:
    return (
        FSL_DERIV
        / f"sub-{sub}"
        / f"ses-{ses}"
        / f"L1_sub-{sub}_ses-{ses}_task-{task}_model-1_type-act_run-{run}_space-mni_{echo}_{cnfd}_unsmoothed_FIR.feat"
    )


def main() -> None:
    if not EV_BASE.exists():
        print(f"ERROR: EVFiles-FIR directory not found:\n  {EV_BASE}")
        return

    skipped_sessions = set()  # track which (sub,ses) we actually encountered + skipped

    missing_inputs = []   # (what, path, context)
    missing_cope1  = []   # cope1 paths (inputs exist => cope1 expected but missing)

    expected_by = defaultdict(int)   # (task, echo, cnfd) -> expected
    present_by  = defaultdict(int)   # (task, echo, cnfd) -> present

    # Walk EVFiles-FIR tree and honor SKIP for any session we encounter.
    for sub_dir in sorted(EV_BASE.glob("sub-*")):
        sub = _parse_id(sub_dir.name, "sub")
        if not sub:
            continue

        for ses_dir in sorted(sub_dir.glob("ses-*")):
            ses = _parse_id(ses_dir.name, "ses")
            if not ses:
                continue

            if (sub, ses) in SKIP:
                skipped_sessions.add((sub, ses))
                continue

            for task in TASKS:
                task_dir = ses_dir / task
                if not task_dir.exists():
                    missing_inputs.append((
                        "EV task directory",
                        task_dir,
                        f"sub-{sub} ses-{ses} task-{task}",
                    ))
                    continue

                for run_dir in sorted(task_dir.glob("run-*")):
                    run = _parse_id(run_dir.name, "run")
                    if not run:
                        continue

                    for echo in ECHOES:
                        bold = _bold_path(sub, ses, task, run, echo)
                        bold_ok = bold.exists()

                        for cnfd in CNFDS:
                            conf = _confounds_path(sub, ses, task, run, cnfd)
                            conf_ok = conf.exists()

                            # Only expect cope1 when prerequisites exist
                            if not bold_ok:
                                missing_inputs.append((
                                    "BOLD",
                                    bold,
                                    f"sub-{sub} ses-{ses} task-{task} run-{run} echo={echo}",
                                ))
                                continue

                            if not conf_ok:
                                missing_inputs.append((
                                    "Confounds TSV",
                                    conf,
                                    f"sub-{sub} ses-{ses} task-{task} run-{run} echo={echo} cnfds={cnfd}",
                                ))
                                continue

                            feat = _feat_dir(sub, ses, task, run, echo, cnfd)
                            cope1 = feat / "stats" / "cope1.nii.gz"

                            expected_by[(task, echo, cnfd)] += 1
                            if cope1.exists():
                                present_by[(task, echo, cnfd)] += 1
                            else:
                                missing_cope1.append(str(cope1))

    # De-dup input-missing lines (they repeat across cnfd in some cases)
    uniq_inputs = {}
    for what, path, ctx in missing_inputs:
        uniq_inputs[(what, str(path), ctx)] = None
    missing_inputs = sorted(uniq_inputs.keys(), key=lambda x: (x[0], x[2], x[1]))
    missing_cope1 = sorted(set(missing_cope1))

    # --------------------------- Print summary ---------------------------
    print("\n=== FIR L1 QC: cope1 presence (mid + sharedreward) ===")
    print(f"Project root: {ROOT_DIR}")
    print(f"EV base:       {EV_BASE}")
    print(f"FSL deriv:     {FSL_DERIV}\n")

    if SKIP:
        # Show which SKIP entries were actually relevant given EVFiles-FIR on disk
        actually_skipped = sorted(skipped_sessions)
        not_seen = sorted(set(SKIP) - skipped_sessions)

        print("Skip list (sub,ses):")
        print("  " + ", ".join([f"{s}:{se}" for s, se in sorted(SKIP)]))
        if actually_skipped:
            print("Skipped (present in EVFiles-FIR and ignored):")
            print("  " + ", ".join([f"{s}:{se}" for s, se in actually_skipped]))
        if not_seen:
            print("Skip entries not found under EVFiles-FIR (nothing to skip):")
            print("  " + ", ".join([f"{s}:{se}" for s, se in not_seen]))
        print("")

    print("Expected -> Present (by task/echo/confounds):")
    if expected_by:
        for key in sorted(expected_by.keys()):
            task, echo, cnfd = key
            exp = expected_by[key]
            pres = present_by.get(key, 0)
            print(f"  {task:12s} {echo:11s} {cnfd:13s}  {exp:4d} -> {pres:4d}   (missing {exp - pres})")
    else:
        print("  (No expected outputs inferred—check EVFiles-FIR tree and inputs.)")
    print("")

    if missing_inputs:
        print("Missing prerequisites (so cope1 is not expected for these rows):")
        for what, path_str, ctx in missing_inputs:
            print(f"  - {what}: {path_str}")
            print(f"      {ctx}")
        print("")
    else:
        print("No missing prerequisites detected.\n")

    if missing_cope1:
        print(f"Missing cope1 (inputs exist; cope1 expected): {len(missing_cope1)}")
        for p in missing_cope1:
            print(f"  - {p}")
        print("")
    else:
        print("No missing cope1 detected.\n")


if __name__ == "__main__":
    main()

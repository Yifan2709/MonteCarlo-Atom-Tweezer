param(
    [string]$Python = 'tmp/long_distance_venv/Scripts/python.exe',
    [int]$Threads = 8
)
$ErrorActionPreference = 'Stop'
$studyPath = $PSScriptRoot
$repoPath = Split-Path (Split-Path $studyPath -Parent) -Parent
Set-Location -LiteralPath $repoPath
$env:PYTHONPATH = "$repoPath/tmp/audit_deps;$repoPath/tmp/continuous_deps;$repoPath/simulation/level2_joint_transfer/src"
$env:NUMBA_NUM_THREADS = "$Threads"
& $Python "$studyPath/validate.py"
if ($LASTEXITCODE -ne 0) { throw 'Physics validation failed' }
& $Python "$studyPath/production_noise_check.py"
if ($LASTEXITCODE -ne 0) { throw 'Production noise validation failed' }
& $Python "$studyPath/papers_and_optics.py"
if ($LASTEXITCODE -ne 0) { throw 'External reference extraction failed' }
& $Python -u "$studyPath/run_study.py" native native_steps coarse mechanisms repeat sensitivity refine confirm convergence operations_confirm repeat_distances repeat_confirm final_single
if ($LASTEXITCODE -ne 0) { throw 'Study run failed' }
& $Python "$studyPath/analyze.py"
if ($LASTEXITCODE -ne 0) { throw 'Analysis failed' }
& $Python -u "$studyPath/run_study.py" boundary_steps edge_robustness handoff_controls repeat_steps
if ($LASTEXITCODE -ne 0) { throw 'Boundary controls failed' }
& $Python "$studyPath/analyze.py"
if ($LASTEXITCODE -ne 0) { throw 'Final analysis failed' }
& $Python "$studyPath/diagnostics.py"
if ($LASTEXITCODE -ne 0) { throw 'Diagnostics failed' }
& $Python "$studyPath/supplement.py"
if ($LASTEXITCODE -ne 0) { throw 'Evidence tables failed' }
& $Python "$studyPath/verify_artifacts.py"
if ($LASTEXITCODE -ne 0) { throw 'Raw-record verification failed' }
& $Python "$studyPath/package_evidence.py"
if ($LASTEXITCODE -ne 0) { throw 'Evidence packaging failed' }

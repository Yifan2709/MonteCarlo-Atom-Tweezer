param(
    [ValidateSet('diagnostics','screen','refine','final','analyze','benchmarks')]
    [string]$Phase='final',
    [string]$Python='',
    [string]$NewRunDirectory=''
)
$ErrorActionPreference='Stop'
$studyDir=$PSScriptRoot
$workspaceDir=(Resolve-Path (Join-Path $studyDir '../..')).Path
if (-not $Python) {
    $Python=Join-Path $studyDir '.venv/Scripts/python.exe'
    if (-not (Test-Path -LiteralPath $Python)) { $Python='python' }
}
Push-Location $workspaceDir
try {
    $queues=switch ($Phase) {
        'diagnostics' { @('diagnostic_queue.json') }
        'screen' { @('scan_queue.json','aod_queue.json','waist_queue.json','lens_queue.json','opt3_queue.json') }
        'refine' { @('refine0_queue.json','refine3_queue.json','selection_queue.json') }
        'final' { @('final_queue.json','validation_checks_queue.json','sensitivity_queue.json') }
        default { @() }
    }
    foreach ($queue in $queues) {
        $queuePath=Join-Path $studyDir $queue
        if (-not (Test-Path -LiteralPath $queuePath)) { continue }
        $invokeArgs=@((Join-Path $studyDir 'run_study.py'),$queuePath)
        if ($NewRunDirectory) { $invokeArgs+=@('--out-dir',$NewRunDirectory) }
        & $Python @invokeArgs
        if ($LASTEXITCODE -ne 0) { throw "Calculation failed: $queue" }
    }
    if ($Phase -eq 'final') {
        $extraArgs=@((Join-Path $studyDir 'extra_timestep.py'))
        if ($NewRunDirectory) { $extraArgs+=@('--runs-dir',$NewRunDirectory) }
        & $Python @extraArgs
        if ($LASTEXITCODE -ne 0) { throw 'Additional timestep check failed' }
    }
    if ($Phase -eq 'analyze') {
        & $Python (Join-Path $studyDir 'diagnostic_plots.py')
        if ($LASTEXITCODE -ne 0) { throw 'Diagnostic plots failed' }
        & $Python (Join-Path $studyDir 'analyze_study.py')
        if ($LASTEXITCODE -ne 0) { throw 'Analysis failed' }
        & $Python (Join-Path $studyDir 'validate_results.py')
        if ($LASTEXITCODE -ne 0) { throw 'Numerical validation analysis failed' }
        & $Python (Join-Path $studyDir 'write_report.py')
        if ($LASTEXITCODE -ne 0) { throw 'Report failed' }
    }
    if ($Phase -eq 'benchmarks') {
        foreach ($script in @('noise_benchmarks.py','gaussian_benchmark.py','waveform_audit.py','static_hold.py')) {
            & $Python (Join-Path $studyDir $script)
            if ($LASTEXITCODE -ne 0) { throw "Benchmark failed: $script" }
        }
    }
} finally { Pop-Location }

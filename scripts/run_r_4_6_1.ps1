[CmdletBinding()]
param(
    [Parameter(Mandatory = $true, Position = 0)]
    [string]$Script,

    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$ScriptArgs
)

$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$projectLibrary = Join-Path $projectRoot '.r-lib\4.6.1'

$rscriptCandidates = @()
if ($env:R461_RSCRIPT) {
    $rscriptCandidates += $env:R461_RSCRIPT
}
$rscriptCommand = Get-Command Rscript.exe -ErrorAction SilentlyContinue
if ($rscriptCommand) {
    $rscriptCandidates += $rscriptCommand.Source
}
if ($env:ProgramFiles) {
    $rscriptCandidates += Join-Path $env:ProgramFiles 'R\R-4.6.1\bin\x64\Rscript.exe'
}
$rscript = $rscriptCandidates |
    Where-Object { $_ -and (Test-Path -LiteralPath $_ -PathType Leaf) } |
    Select-Object -First 1
if (-not $rscript) {
    throw 'R 4.6.1 Rscript was not found. Set R461_RSCRIPT or add Rscript.exe to PATH.'
}

$scriptPath = if ([System.IO.Path]::IsPathRooted($Script)) {
    $Script
} else {
    Join-Path $projectRoot $Script
}
if (-not (Test-Path -LiteralPath $scriptPath -PathType Leaf)) {
    throw "R script was not found: $scriptPath"
}

# Git Bash and some inherited shells set C.UTF-8, which is not a valid Windows
# locale name for R. Remove only in this launcher process and its R child.
foreach ($name in @(
    'LC_ALL', 'LC_COLLATE', 'LC_CTYPE', 'LC_MONETARY',
    'LC_NUMERIC', 'LC_TIME', 'LANG', 'LANGUAGE'
)) {
    Remove-Item -Path "Env:$name" -ErrorAction SilentlyContinue
}

# Keep this pipeline isolated from the per-user R library. In particular, the
# latter may contain Windows arrow builds that can crash Rscript during process
# teardown. Revision scripts use the pinned project-local nanoparquet library.
$env:R_LIBS_USER = $projectLibrary

$exitCode = 1
Push-Location -LiteralPath $projectRoot
try {
    & $rscript --vanilla $scriptPath @ScriptArgs
    $exitCode = $LASTEXITCODE
} finally {
    Pop-Location
}
exit $exitCode

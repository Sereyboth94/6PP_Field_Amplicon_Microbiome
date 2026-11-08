#Requires -Version 5.1
<# Copy exactly twelve processed input tables; original folders remain unchanged.
   Default SourceRoot is the parent of this repository, matching the user's layout.
   -CheckOnly reads sources and hashes without copying anything.
#>
[CmdletBinding()]
param(
    [string]$SourceRoot = (Split-Path -Parent (Split-Path -Parent $PSScriptRoot)),
    [switch]$CheckOnly
)
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot 'Repo_Common.ps1')
$repoRoot = (Get-Item -LiteralPath (Split-Path -Parent $PSScriptRoot)).FullName
$sourceAbsolute = (Get-Item -LiteralPath $SourceRoot).FullName
if ($sourceAbsolute.TrimEnd('\', '/') -eq $repoRoot.TrimEnd('\', '/')) {
    throw 'SourceRoot must be the original github folder outside this repository.'
}
$relativePaths = @(Get-AmpliconInputPaths)
$missing = @($relativePaths | Where-Object {
    -not (Test-Path -LiteralPath (Join-Path $sourceAbsolute $_) -PathType Leaf)
})
if ($missing.Count -gt 0) { throw ('Original input files are missing. No files were copied:' + [Environment]::NewLine + ($missing -join [Environment]::NewLine)) }
$checks = foreach ($relative in $relativePaths) {
    $source = Join-Path $sourceAbsolute $relative
    $item = Get-Item -LiteralPath $source
    if ($item.Length -eq 0) { throw "Original input is empty: $source" }
    [PSCustomObject]@{ File = $relative; Bytes = $item.Length; SHA256 = (Get-FileHash -LiteralPath $source -Algorithm SHA256).Hash }
}
if ($CheckOnly) {
    $checks | Format-Table File, Bytes -AutoSize
    Write-Host 'All twelve original input files are readable. CheckOnly made no changes.'
    return
}
$stage = Join-Path $repoRoot ('.input_staging_' + [guid]::NewGuid().ToString('N'))
$backupRoot = Join-Path $repoRoot ('.input_backups/' + (Get-Date -Format 'yyyyMMdd_HHmmss') + '_' + [guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Path $stage | Out-Null
try {
    # Verify every staged copy before changing any existing repository input.
    foreach ($check in $checks) {
        $staged = Join-Path $stage $check.File
        New-Item -ItemType Directory -Path (Split-Path -Parent $staged) -Force | Out-Null
        Copy-Item -LiteralPath (Join-Path $sourceAbsolute $check.File) -Destination $staged
        if ((Get-FileHash -LiteralPath $staged -Algorithm SHA256).Hash -ne $check.SHA256) { throw "Copied checksum mismatch: $($check.File)" }
    }
    foreach ($check in $checks) {
        $target = Join-Path $repoRoot $check.File
        New-Item -ItemType Directory -Path (Split-Path -Parent $target) -Force | Out-Null
        if (Test-Path -LiteralPath $target -PathType Leaf) {
            $oldHash = (Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash
            if ($oldHash -eq $check.SHA256) { continue }
            $backup = Join-Path $backupRoot $check.File
            New-Item -ItemType Directory -Path (Split-Path -Parent $backup) -Force | Out-Null
            Copy-Item -LiteralPath $target -Destination $backup
            if ((Get-FileHash -LiteralPath $backup -Algorithm SHA256).Hash -ne $oldHash) { throw "Backup verification failed: $target" }
        }
        Copy-Item -LiteralPath (Join-Path $stage $check.File) -Destination $target -Force
        if ((Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash -ne $check.SHA256) { throw "Installed checksum mismatch: $target" }
    }
    $csv = [string[]]@($checks | ConvertTo-Csv -NoTypeInformation)
    [IO.File]::WriteAllLines((Join-Path $repoRoot 'data_checksums.csv'), $csv, [Text.UTF8Encoding]::new($false))
    Write-Host 'Imported and SHA256-verified all twelve processed tables.' -ForegroundColor Green
    Write-Host 'Original folders were not changed. Input checksums: data_checksums.csv'
    if (Test-Path -LiteralPath $backupRoot) { Write-Host ('Replaced repository input backups: ' + $backupRoot) }
} finally {
    if (Test-Path -LiteralPath $stage) { Remove-Item -LiteralPath $stage -Recurse -Force }
}

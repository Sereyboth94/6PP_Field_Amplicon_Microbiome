# Shared helpers used by the local import and GitHub publication scripts.
function Get-AmpliconInputPaths {
    $groups = @(
        @('Pasture', 'Bacteria', 'Metadat_Bacteria_Field.txt'),
        @('Pasture', 'Fungi', 'Metadat_Fungi_Pasture.txt'),
        @('Wheat', 'Bacteria', 'Metadat_Bacteria_Field.txt'),
        @('Wheat', 'Fungi', 'Metadat_Fungal_Wheat.txt')
    )
    foreach ($group in $groups) {
        foreach ($name in @('ASV_filtered.txt', 'Taxonomy_filtered.txt', $group[2])) {
            $group[0] + '/' + $group[1] + '/' + $name
        }
    }
}

function Find-AmpliconApplication {
    param([string]$Name, [string[]]$RelativePaths)
    $command = Get-Command $Name -CommandType Application -ErrorAction SilentlyContinue | Select-Object -First 1
    if ($command) { return $command.Source }
    foreach ($base in @($env:ProgramFiles, ${env:ProgramFiles(x86)}, $env:LOCALAPPDATA)) {
        if ([string]::IsNullOrWhiteSpace($base)) { continue }
        foreach ($relative in $RelativePaths) {
            $candidate = Join-Path $base $relative
            if (Test-Path -LiteralPath $candidate -PathType Leaf) { return $candidate }
        }
    }
    return $null
}

function Invoke-AmpliconNative {
    param([string]$Executable, [string[]]$ArgumentList)
    & $Executable @ArgumentList
    if ($LASTEXITCODE -ne 0) { throw ('Command failed with exit code {0}: {1}' -f $LASTEXITCODE, [IO.Path]::GetFileName($Executable)) }
}

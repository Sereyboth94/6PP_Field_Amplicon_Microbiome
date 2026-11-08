#Requires -Version 5.1
<# Create and push the prepared repository to the authenticated GitHub account.
   Default: Sereyboth94/6PP_Field_Amplicon_Microbiome, private.
   A pre-existing remote repository is left unchanged; no force push is used.
#>
[CmdletBinding()]
param(
    [ValidatePattern('^[A-Za-z0-9-]+$')][string]$Owner = 'Sereyboth94',
    [ValidatePattern('^[A-Za-z0-9_.-]+$')][string]$RepositoryName = '6PP_Field_Amplicon_Microbiome',
    [ValidateSet('Private', 'Public')][string]$Visibility = 'Private'
)
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot 'Repo_Common.ps1')
$repoRoot = (Get-Item -LiteralPath (Split-Path -Parent $PSScriptRoot)).FullName
$inputPaths = @(Get-AmpliconInputPaths)
foreach ($relative in $inputPaths) {
    $path = Join-Path $repoRoot $relative
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) { throw "Missing input: $relative. Run Import_Original_Data.ps1 before publishing." }
    if ((Get-Item -LiteralPath $path).Length -eq 0) { throw "Empty input: $relative" }
}
$manifestPath = Join-Path $repoRoot 'data_checksums.csv'
if (-not (Test-Path -LiteralPath $manifestPath -PathType Leaf)) { throw 'Run Import_Original_Data.ps1 to create the verified input manifest first.' }
$manifest = @(Import-Csv -LiteralPath $manifestPath)
if ($manifest.Count -ne 12) { throw 'The input manifest must contain exactly twelve files.' }
foreach ($relative in $inputPaths) {
    $entry = @($manifest | Where-Object { $_.File -eq $relative })
    if ($entry.Count -ne 1 -or (Get-FileHash -LiteralPath (Join-Path $repoRoot $relative) -Algorithm SHA256).Hash -ne $entry[0].SHA256) {
        throw "Input changed after import: $relative. Review it and rerun the import before publishing."
    }
}
$git = Find-AmpliconApplication 'git' @('Git/cmd/git.exe', 'Programs/Git/cmd/git.exe')
$gh = Find-AmpliconApplication 'gh' @('GitHub CLI/gh.exe', 'Programs/GitHub CLI/gh.exe')
if (-not $git) { throw 'Install Git for Windows (https://git-scm.com/download/win), then rerun this script.' }
if (-not $gh) { throw 'Install GitHub CLI (https://cli.github.com/), then rerun this script. Your prepared files are retained.' }
& $gh auth status --hostname github.com
if ($LASTEXITCODE -ne 0) {
    Invoke-AmpliconNative $gh @('auth', 'login', '--hostname', 'github.com', '--git-protocol', 'https', '--web', '--scopes', 'workflow')
}
$login = (Invoke-AmpliconNative $gh @('api', 'user', '--jq', '.login') | Out-String).Trim()
if ($login -ine $Owner) { throw "Authenticated as $login, but the requested owner is $Owner. Switch your active GitHub CLI account and rerun." }
$fullName = $Owner + '/' + $RepositoryName
$oldPreference = $ErrorActionPreference
$ErrorActionPreference = 'Continue'
try {
    $existing = & $gh repo view $fullName --json nameWithOwner --jq '.nameWithOwner' 2>$null
    $viewExit = $LASTEXITCODE
} finally { $ErrorActionPreference = $oldPreference }
if ($viewExit -eq 0) { throw "Repository $fullName already exists. It was not changed. Use another -RepositoryName for a new repository." }
if (-not (Test-Path -LiteralPath (Join-Path $repoRoot '.git'))) {
    Invoke-AmpliconNative $git @('-C', $repoRoot, 'init', '-b', 'main')
}
$remotes = @(& $git -C $repoRoot remote)
if ($remotes.Count -gt 0) { throw 'This local repository already has a remote. Review that remote before creating another hosted repository.' }
$gitName = & $git -C $repoRoot config user.name
$gitEmail = & $git -C $repoRoot config user.email
if (-not $gitName) { Invoke-AmpliconNative $git @('-C', $repoRoot, 'config', '--local', 'user.name', $login) }
if (-not $gitEmail) {
    $numericId = (Invoke-AmpliconNative $gh @('api', 'user', '--jq', '.id') | Out-String).Trim()
    if ($numericId -notmatch '^[0-9]+$') { throw 'Could not determine the authenticated GitHub user ID.' }
    Invoke-AmpliconNative $git @('-C', $repoRoot, 'config', '--local', 'user.email', ($numericId + '+' + $login + '@users.noreply.github.com'))
}
# Explicit staging excludes local outputs, raw reads, and unrelated source files.
$pathsToStage = @('README.md', '.gitignore', '.gitattributes', 'CITATION.cff',
    '6PP_Field_Amplicon_Microbiome.Rproj', 'Microbiome_Figures_Scripts_and_Tables',
    'R', 'tools', 'tests', 'docs', '.github', 'data_checksums.csv') + $inputPaths
Invoke-AmpliconNative $git (@('-C', $repoRoot, 'add', '--') + $pathsToStage)
$staged = @(& $git -C $repoRoot diff --cached --name-only)
if ($staged.Count -gt 0) {
    Invoke-AmpliconNative $git @('-C', $repoRoot, 'commit', '-m', 'Add reproducible pasture and wheat microbiome figure workflows')
}
Invoke-AmpliconNative $git @('-C', $repoRoot, 'rev-parse', '--verify', 'HEAD') | Out-Null
$visibilityFlag = '--' + $Visibility.ToLowerInvariant()
Invoke-AmpliconNative $gh @('repo', 'create', $fullName, $visibilityFlag, '--source', $repoRoot, '--remote', 'origin', '--push',
    '--description', 'Processed-input R workflows for pasture and wheat microbiome diversity, composition, and manuscript figures')
$url = (Invoke-AmpliconNative $gh @('repo', 'view', $fullName, '--json', 'url', '--jq', '.url') | Out-String).Trim()
Write-Host ('Created and pushed ' + $Visibility.ToLowerInvariant() + ' repository: ' + $url) -ForegroundColor Green

# Dot-source this file: . ./tools/use_gcc.ps1 -GccBin 'D:/.../mingw64/bin'
# Applies to the current PowerShell session only.
param([string]$GccBin = 'D:/Code_Configuration/gcc-16.2.0/mingw64/bin')
$ErrorActionPreference = 'Stop'
$compilerDirectory = (Resolve-Path -LiteralPath $GccBin).Path
foreach ($compiler in @('gcc.exe', 'g++.exe', 'mingw32-make.exe')) {
    if (-not (Test-Path -LiteralPath (Join-Path $compilerDirectory $compiler))) {
        throw "Missing compiler tool: $compiler"
    }
}
$env:PATH = $compilerDirectory + ';' + $env:PATH
$env:CC = Join-Path $compilerDirectory 'gcc.exe'
$env:CXX = Join-Path $compilerDirectory 'g++.exe'
$fortranCompiler = Join-Path $compilerDirectory 'gfortran.exe'
if (Test-Path -LiteralPath $fortranCompiler) { $env:FC = $fortranCompiler }
$env:AEROACOUSTICS_RUNTIME_DIRS = $compilerDirectory
& $env:CC --version | Select-Object -First 1
& $env:CXX --version | Select-Object -First 1
if (Test-Path -LiteralPath $fortranCompiler) { & $fortranCompiler --version | Select-Object -First 1 }

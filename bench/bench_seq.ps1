<#
.SYNOPSIS
  Mjerenje sekvencijalne verzije na svim veličinama mreže.

.DESCRIPTION
  Pokreće erosion_cli za 256², 512², 1024² i 2048² (puni korak simulacije:
  voda, erozija, transport, isparavanje), svaku veličinu više puta, vezano za
  jedno jezgro (zadano CPU 0 = P-jezgro na hibridnim Intel procesorima).
  Sirovi rezultati idu u CSV, a na kraju se ispisuje medijan po veličini.

.EXAMPLE
  powershell -ExecutionPolicy Bypass -File bench\bench_seq.ps1
  powershell -ExecutionPolicy Bypass -File bench\bench_seq.ps1 -Repeats 5 -Exe build\erosion_cli.exe
#>
param(
    [string]$Exe       = "build\erosion_cli.exe",
    [string]$OutCsv    = "bench\results\seq.csv",
    [int]   $Repeats   = 3,
    [long]  $Affinity  = 0x1,   # maska jezgara; 0x1 = CPU 0 (P-jezgro na i7-1355U)
    [switch]$AllowBattery
)

$ErrorActionPreference = "Stop"

if (-not (Test-Path $Exe)) {
    throw "Ne postoji $Exe. Prvo: cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release; cmake --build build"
}

# Na bateriji Windows usporava procesor (izmjereno ~2x sporije, nestabilno).
$battery = Get-CimInstance -Namespace root/wmi -ClassName BatteryStatus -ErrorAction SilentlyContinue |
           Select-Object -First 1
if ($battery -and -not $battery.PowerOnline -and -not $AllowBattery) {
    throw "Laptop je na bateriji: rezultati ne bi bili pouzdani. Ukljucite punjac (ili -AllowBattery za probu)."
}

$sizes = @(
    @{ Size = 256;  Steps = 1000 },
    @{ Size = 512;  Steps = 300  },
    @{ Size = 1024; Steps = 80   },
    @{ Size = 2048; Steps = 20   }
)

New-Item -ItemType Directory -Force (Split-Path $OutCsv) | Out-Null
$tmpDir = Join-Path ([IO.Path]::GetTempPath()) "prs_bench"
New-Item -ItemType Directory -Force $tmpDir | Out-Null
$exePath = (Resolve-Path $Exe).Path
$csvPath = if ([IO.Path]::IsPathRooted($OutCsv)) { $OutCsv } else { Join-Path (Get-Location) $OutCsv }
$csvPath = [IO.Path]::GetFullPath($csvPath)

foreach ($cfg in $sizes) {
    foreach ($rep in 1..$Repeats) {
        $cliArgs = "--size $($cfg.Size) --scenario layer --depth 1 --rain 0.02 --steps $($cfg.Steps) " +
                   "--out terrain.pgm --water-out water.pgm --csv `"$csvPath`" --tag seq"
        $p = Start-Process -FilePath $exePath -ArgumentList $cliArgs -WorkingDirectory $tmpDir `
                           -NoNewWindow -PassThru -RedirectStandardOutput (Join-Path $tmpDir "out.txt")
        $null = $p.Handle  # bez ovoga je $p.ExitCode prazan nakon WaitForExit (poznata zamka PowerShella)
        try { $p.ProcessorAffinity = [IntPtr]$Affinity } catch { Write-Warning "Ne mogu postaviti affinity: $_" }
        $p.WaitForExit()
        if ($p.ExitCode -ne 0) { throw "erosion_cli nije uspio (kod $($p.ExitCode)) za $($cfg.Size)^2" }
        Write-Host ("{0,5}^2  ponavljanje {1}/{2} gotovo" -f $cfg.Size, $rep, $Repeats)
    }
}

# Medijan po veličini, samo za ovo pokretanje (zadnjih N redova po veličini).
$rows = Import-Csv $csvPath | Where-Object { $_.tag -eq "seq" }
function Median([double[]]$v) { $s = $v | Sort-Object; $n = $s.Count; if ($n % 2) { $s[[int][math]::Floor($n / 2)] } else { ($s[$n / 2 - 1] + $s[$n / 2]) / 2 } }

Write-Host ""
Write-Host "Medijan (ms/korak) po fazama:"
"{0,6} {1,10} {2,9} {3,8} {4,9} {5,8} {6,8} {7,8} {8,8}" -f "N", "ukupno", "Mcel/s", "voda+", "fluks", "transp", "voda", "erozija", "ispar."
foreach ($cfg in $sizes) {
    $r = @($rows | Where-Object { [int]$_.width -eq $cfg.Size } | Select-Object -Last $Repeats)
    if ($r.Count -eq 0) { continue }
    "{0,6} {1,10:N3} {2,9:N1} {3,8:N3} {4,9:N3} {5,8:N3} {6,8:N3} {7,8:N3} {8,8:N3}" -f `
        $cfg.Size,
        (Median ($r | ForEach-Object { [double]$_.ms_per_step })),
        (Median ($r | ForEach-Object { [double]$_.mcells_per_s })),
        (Median ($r | ForEach-Object { [double]$_.add_water })),
        (Median ($r | ForEach-Object { [double]$_.flux })),
        (Median ($r | ForEach-Object { [double]$_.transport })),
        (Median ($r | ForEach-Object { [double]$_.water })),
        (Median ($r | ForEach-Object { [double]$_.erosion })),
        (Median ($r | ForEach-Object { [double]$_.evaporate }))
}
Write-Host ""
Write-Host "Sirovi rezultati: $OutCsv"

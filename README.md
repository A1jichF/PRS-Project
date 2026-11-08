# Paralelna simulacija plitke vode i erozije terena (P11)

Projekat iz predmeta **Paralelni računarski sistemi**, tema **P11: Shallow-water and Terrain Erosion**.

Cilj je implementirati simulaciju hidrauličke erozije terena zasnovanu na jednačinama plitke vode (shallow-water equations), i to po **virtual pipe** modelu (Mei et al. 2007). Simulacija radi na 2D mreži visina (heightmap). Prvo pravimo sekvencijalnu verziju, a zatim je paralelizujemo na CPU-u pomoću OpenMP-a i SIMD-a (AVX2).

## Dokumentacija

Folder [Izvori/](Izvori/) sadrži naučne radove, linkove i ostale dokumente projekta.

- **[Pregled metoda erozije/](Izvori/Pregled%20metoda%20erozije/)**: pregled oblasti i drugih pristupa simulaciji erozije, da bismo opravdali zašto ih ne biramo.
  - *A Review of Digital Terrain Modeling* (Galin et al., Eurographics STAR 2019): opšti pregled modeliranja terena, uključujući eroziju.
  - *Hydraulic Erosion* (Beneš et al., CAVW 2006): 3D (voxel/Navier–Stokes) pristup eroziji.
  - *Hydraulic Erosion Using Smoothed Particle Hydrodynamics* (Krištof et al., Eurographics 2009): čestični (SPH) pristup.
- **[Zasto Shallow-Water/](Izvori/Zasto%20Shallow-Water/)**: radovi koji su osnova za izabrani model.
  - *Fast Hydraulic Erosion Simulation and Visualization on GPU* (Mei, Decaudin, Hu, PG 2007): **osnova algoritma** (virtual pipe model + erozija/depozicija).
  - *Interactive Terrain Modeling Using Hydraulic Erosion* (Šťava et al., SCA 2008): proširenje na više slojeva materijala, odronjavanje obala i podjelu na blokove (tiling).
  - *A Second-Order Well-Balanced Positivity Preserving Central-Upwind Scheme for the Saint-Venant System* (Kurganov, Petrova 2007): numerički tačnija shema za shallow-water jednačine, kao poređenje s pipe modelom.
- **[Paralelizacije/](Izvori/Paralelizacije/)**: paralelne (GPU/CPU) implementacije shallow-water i erozije.
  - *Simulation and Visualization of the Saint-Venant System using GPUs* (Brodtkorb et al., 2010)
  - *Fast Hydraulic and Thermal Erosion on the GPU* (Jákó, CESCG 2011)
  - *3D Real-Time Hydraulic Erosion Simulation using Multi-Layered Heightmaps* (Nilles et al., VMV 2024)
  - Šira referenca (PDF ručno): *Efficient Shallow Water Simulations on GPUs* (Brodtkorb, Sætra, Altinakar 2012) i *TRITON: A Multi-GPU Open Source 2D Hydrodynamic Flood Model* (Morales-Hernández et al. 2021).
- **[linkovi-implementacije.md](Izvori/linkovi-implementacije.md)**: GitHub repozitoriji sa implementacijama i ostali korisni linkovi.

## Početni plan implementacije

C++20 i CMake. Teren je mreža `N × N` ćelija, a svaka ćelija čuva: visinu terena `b`, visinu vode `d`, sediment `s`, fluks prema 4 susjeda `f_L, f_R, f_T, f_B` i brzinu `(u, v)`. Podaci su u **SoA** rasporedu (posebni nizovi po veličini) zbog SIMD-a.

1. **Sekvencijalna verzija.** Jedan korak simulacije (po Mei et al.):
   1. dodavanje vode (kiša / izvor),
   2. računanje fluksa kroz virtualne cijevi,
   3. ažuriranje visine vode i polja brzine,
   4. erozija i depozicija (kapacitet nošenja sedimenta),
   5. transport sedimenta (semi-Lagrangian advekcija; zadano konzervativni transport kroz cijevi, vidi ispod),
   6. isparavanje,
   7. (opcionalno) termalna erozija.

   Ulaz: generisani teren (Perlin/fBm šum) ili DEM. Izlaz: heightmap kao PGM/PNG.
   Metrike: vrijeme po koraku, broj ažuriranih ćelija u sekundi (cells/s), očuvanje mase vode i materijala.
2. **SIMD.** Kernele fluksa i erozije pišemo u tri verzije: skalarna, `#pragma omp simd` i ručni AVX2/FMA intrinsics.
3. **OpenMP.** Paralelizacija svakog koraka po redovima mreže (`parallel for`, poređenje `static` i `dynamic` rasporeda). Koristimo double-buffering da nema data race-a između čitanja susjeda i pisanja. Advekcija sedimenta je nepravilan pristup memoriji, pa nju posebno mjerimo. Probat ćemo i blokovsku (tiled) podjelu zbog keša.
4. **Mjerenja.** Veličine mreže ∈ {256², 512², 1024², 2048²}, niti ∈ {1, 2, 4, 8}, ubrzanje i efikasnost u odnosu na sekvencijalnu verziju, te strong i weak scaling.

## Build i pokretanje

Potrebni su CMake ≥ 3.20 i C++20 kompajler (MSVC, GCC ili Clang). Na Windowsu u Visual Studio Installeru uključite workload **Desktop development with C++** (sadrži i CMake i Ninja).

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure

# generiše fBm teren 512 x 512 i upisuje 16-bitni PGM
build/erosion_cli --size 512 --seed 42 --out terrain.pgm
build/erosion_cli --help

# sloj vode 1 m preko fBm terena, 2000 koraka -> voda se skupi u dolinama
build/erosion_cli --size 256 --scenario layer --depth 1 --steps 2000 --water-out water.pgm

# dam-break sa snimkom vode svakih 150 koraka (db_000000.pgm, db_000150.pgm, ...)
build/erosion_cli --size 256 --cell-size 1 --scenario dambreak --steps 600 --snapshot-every 150 --snapshot-prefix db

# kiša i izvor vode (x, y, poluprečnik, količina u sekundi)
build/erosion_cli --size 256 --rain 0.01 --source 128 40 4 5 --steps 3000

# erozija pod kišom: mapa promjene terena (tamno = erozija, svijetlo = taloženje) i sediment
build/erosion_cli --size 256 --rain 0.02 --steps 5000 --diff-out diff.pgm --sediment-out sediment.pgm

# ista simulacija bez erozije, ili s drugim konstantama
build/erosion_cli --size 256 --rain 0.02 --steps 5000 --no-erosion
build/erosion_cli --size 256 --rain 0.02 --steps 5000 --kc 0.3 --ks 0.3 --kd 1 --kdmax 0.5
```

Na Windowsu s MSYS2 (UCRT64) GCC-om programi se linkuju statički (`-static`), pa rade i kad MSYS2 nije na PATH-u.

### Mjerenje vremena na hibridnim procesorima

Intel procesori od 12. generacije imaju jaka (P) i slabija (E) jezgra. Na i7-1355U je isti sekvencijalni kod na E-jezgru **2,1× sporiji** (82 naspram 39 Mćelija/s na 512²), a Windows sam bira gdje će program raditi. Zato se sekvencijalna mjerenja rade **vezano za P-jezgro**:

```bat
start /affinity 1 /wait /b build\erosion_cli.exe --size 512 --scenario layer --steps 300
```

Svako mjerenje ponoviti barem 3 puta i uzeti medijan, s laptopom **na punjaču** i u režimu "Best performance". Na bateriji je isti kod bio oko 2× sporiji, a vremena su skakala i do 60% između pokretanja.

### Struktura koda

| Putanja | Sadržaj |
|---|---|
| [include/erosion/Grid.h](include/erosion/Grid.h) | Stanje mreže (SoA): `b`, `d`, `s`, `fL/fR/fT/fB`, `u`, `v`, poravnato na 64 B. Indeks `y * width + x`. |
| [include/erosion/Params.h](include/erosion/Params.h) | Parametri simulacije (`dt`, `g`, `A`, `Kc`, `Ks`, `Kd`, `Ke`...). |
| [include/erosion/TerrainGen.h](include/erosion/TerrainGen.h) | Deterministički fBm generator terena. |
| [include/erosion/ImageIO.h](include/erosion/ImageIO.h) | Čitanje i pisanje PGM-a (8 i 16 bita). |
| [include/erosion/Diagnostics.h](include/erosion/Diagnostics.h) | Ukupna masa vode i materijala, min/max, poređenje polja, CFL. |
| [include/erosion/Kernels.h](include/erosion/Kernels.h) | Kerneli po fazama (`seq::addWater`, `seq::computeFlux`, `seq::updateWaterAndVelocity`). |
| [include/erosion/Simulation.h](include/erosion/Simulation.h) | Poziva kernele redom i mjeri vrijeme po fazi. |
| [include/erosion/Scenarios.h](include/erosion/Scenarios.h) | Početna stanja: mirno jezero, punjenje do nivoa, dam-break. |
| [apps/erosion_cli/](apps/erosion_cli/) | Komandna linija: generisanje i učitavanje terena (kasnije i simulacija). |
| [tests/](tests/) | Testovi bez vanjskih zavisnosti. |

Biblioteka `erosion` ne zavisi ni od kakve grafike, pa se kasnije može uključiti u [opengl-forward-renderer](https://github.com/bakir004/opengl-forward-renderer) kao zasebna scena.

## Odstupanja od Mei et al. 2007

| Šta | Mei et al. | Kod nas | Zašto |
|---|---|---|---|
| `Ks`, `Kd` | udio po koraku | brzina [1/s] × `dt` | rezultat ne zavisi od izbora `dt` |
| Kapacitet `C` | `Kc · sin α · \|v\|` | `× lmax(d)` (Šťava et al. 2008) | na ćelijama koje se isprazne u jednom koraku brzina je umjetno `l/dt` |
| Transport sedimenta | semi-Lagrangian `s(x − v·Δt)` | **konzervativno kroz virtualne cijevi** (`--transport flux`, zadano); Mei-jeva metoda ostaje kao `--transport sl` | vidi ispod |

**Transport sedimenta.** Semi-Lagrangian advekcija rješava `∂s/∂t + v·∇s = 0`, dakle tretira vodu kao nestišljivu. U plitkoj vodi to ne vrijedi: tamo gdje se tok usporava i skuplja (dno padine, ulaz u jezero, `∇·v < 0`) sediment nestaje, a to su upravo mjesta gdje bi se trebao taložiti. Na fBm terenu 256² s kišom, nakon 5000 koraka:

| Metoda | Promjena ukupnog materijala | Promjena terena (min / max) |
|---|---|---|
| semi-Lagrangian (Mei) | **−19,8 %** | −70,3 / +31,6 |
| konzervativni kroz cijevi | **+1,2·10⁻⁷** | −29,4 / +21,2 |

Konzervativni transport šalje susjedu isti udio sedimenta kao i vode (isti fluks iz faze 2), pa se masa čuva do greške zaokruživanja, a sediment se taloži nizvodno. Kao i update vode, to je stencil s 4 susjeda (gather, bez konflikta pisanja), pa se jednako dobro paralelizuje. Semi-Lagrangian ostaje u kodu za poređenje, a njegov nepravilan pristup memoriji je zanimljiv za SIMD/OpenMP analizu.

## Način rada

Svaki korak se provjerava prije sljedećeg:
- **Mirno jezero** (ravan teren s vodom) mora ostati mirno, tj. brzina ≈ 0.
- **Ukupna masa vode** (bez kiše i isparavanja) mora ostati konstantna, do greške zaokruživanja.
- **Ukupna masa materijala** (teren + sediment) mora ostati konstantna.
- Visina vode nikad ne smije postati negativna.
- **Dam-break test** (pucanje brane): kvalitativno se poredi sa očekivanim talasom.
- SIMD i OpenMP verzije moraju dati isti rezultat kao sekvencijalna (do tolerancije `1e-5`).

Sve razvijamo na maloj mreži (256²), a konačna mjerenja radimo na velikim.

## Rječnik pojmova

- **Shallow-water jednačine (Saint-Venant)**: pojednostavljene jednačine fluida za slučaj kad je dubina vode mnogo manja od horizontalnih dimenzija.
- **Virtual pipe model**: svaka ćelija je povezana "cijevima" sa 4 susjeda, a fluks kroz cijev ubrzava razlika visina (voda + teren).
- **Heightmap**: 2D mreža visina terena.
- **Fluks**: količina vode koja u jednom koraku pređe iz ćelije u susjednu.
- **Kapacitet nošenja sedimenta (C)**: koliko sedimenta voda može nositi. Ovisi o brzini i nagibu: ako je sediment manji od C, dolazi do erozije, a ako je veći, do depozicije.
- **Semi-Lagrangian advekcija**: sediment u ćeliji se uzima iz tačke `x − v·Δt`, uz bilinearnu interpolaciju.
- **Termalna erozija**: odronjavanje materijala kad je nagib veći od ugla talusa.
- **CFL uslov**: ograničenje veličine vremenskog koraka da bi simulacija bila stabilna.
- **Double-buffering**: dva bafera (staro/novo stanje), tako da paralelne niti čitaju samo staro stanje.
- **Stencil**: obrazac računanja gdje nova vrijednost ćelije ovisi o njoj i njenim susjedima.
- **SoA / AoS**: Structure of Arrays / Array of Structures, raspored podataka u memoriji.
- **OpenMP, SIMD, AVX2, FMA**: paralelizacija po nitima, vektorske instrukcije, 256-bitni skup instrukcija, fused multiply-add.
- **Speedup / efikasnost**: `T_seq / T_par` i `speedup / broj_niti`.

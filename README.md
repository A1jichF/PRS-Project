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

Intel procesori od 12. generacije imaju jaka (P) i slabija (E) jezgra. Na i7-1355U je isti sekvencijalni kod na E-jezgru **2,1× sporiji** (82 naspram 39 Mćelija/s na 512²), a Windows sam bira gdje će program raditi. Zato se sekvencijalna mjerenja rade **vezano za P-jezgro**. To radi skripta [bench/bench_seq.ps1](bench/bench_seq.ps1): pokreće sve 4 veličine mreže po 3 puta na CPU 0, dopisuje rezultate u `bench/results/seq.csv` i ispisuje medijan po fazama:

```powershell
powershell -ExecutionPolicy Bypass -File bench\bench_seq.ps1
```

Skripta odbija da radi kad je laptop na bateriji. Za ručno pokretanje: `start /affinity 1 /wait /b build\erosion_cli.exe ... --csv rezultati.csv`.

Mjerenja raditi s laptopom **na punjaču** i u režimu "Best performance".

### Rezultati: sekvencijalna verzija (`T_seq`)

i7-1355U, CPU 0 (P-jezgro), na punjaču, Windows režim **"Best performance"**, GCC 16.2 `-O2`, scenarij `layer` + kiša 0,02 (sve faze aktivne), medijan 5 ponavljanja. Sirovi podaci: [bench/results/seq.csv](bench/results/seq.csv).

| N | ms/korak | Mćelija/s | fluks | transport | voda + brzina | erozija | dodavanje + isparavanje | rasipanje |
|---|---|---|---|---|---|---|---|---|
| 256² | 0,94 | 69,8 | 0,19 | 0,16 | 0,15 | 0,44 | 0,01 | 2,8 % |
| 512² | 3,54 | 74,1 | 0,76 | 0,69 | 0,61 | 1,43 | 0,06 | 2,6 % |
| 1024² | 12,9 | 81,5 | 3,04 | 3,57 | 2,55 | 3,36 | 0,45 | 3,7 % |
| 2048² | 54,3 | 77,3 | 12,3 | 14,8 | 10,6 | 13,1 | 3,41 | 2,1 % |

U režimu "Balanced" je isto mjerenje bilo 8–11 % sporije i s rasipanjem do 7,3 % ([bench/results/seq_balanced.csv](bench/results/seq_balanced.csv)). Sva kasnija mjerenja (OpenMP, SIMD) rade se u režimu "Best performance".

### Zamka: 4K aliasing

Prvo mjerenje je imalo rasipanje do **68 %**, i to samo u fazi "voda + brzina". Ista faza je, uz bit-identičan ulaz i rezultat, trajala 0,16, 0,38, 0,60 ili 0,82 ms (256²), uvijek isto unutar jednog procesa. Postupnim eliminisanjem (jezgro, ASLR, fizička memorija, vrijeme od starta) pokazalo se da zavisi samo od **međusobnog položaja nizova u memoriji**: zamjena niza `u` svježom kopijom je sporu fazu činila brzom.

Uzrok: sva polja mreže su iste veličine i alocirana jedno za drugim, pa su počinjala na **istom pomaku unutar 4 KB stranice** (`u ≡ fR`, `v ≡ fT`, `fL ≡ fB` po modulu 4096). Petlja upiše `u[i]`, pa u sljedećoj iteraciji čita `fR[i]`. Procesor upoređuje samo nižih 12 bitova adrese, pogrešno zaključi da čitanje zavisi od upisa i čeka. Četiri nivoa brzine (`0,16 + k · 0,22` ms) odgovaraju broju pogođenih izlaznih nizova (`d`, `u`, `v`).

Popravka je u [AlignedAllocator.h](include/erosion/AlignedAllocator.h): svako polje dobija vlastiti pomak unutar stranice (korak 256 B; polja mreže slotovi 0–8, pomoćni baferi simulacije 9–11). Ovo je važno i za OpenMP/SIMD mjerenja, jer bi inače nasumično usporavanje izgledalo kao loše skaliranje. Na bateriji je isti kod bio oko 2× sporiji, a vremena su skakala i do 60% između pokretanja.

### 3D prikaz u opengl-forward-renderer

Simulacija se može gledati uživo u 3D, u forku engina [KerimSikalo/opengl-forward-renderer](https://github.com/KerimSikalo/opengl-forward-renderer/tree/feature/erosion-scene) (grana `feature/erosion-scene`, scena **"Erosion (PRS)"**). Fork treba biti kloniran **pored** ovog foldera:

```
Repozitorij\
  PRS-Project\                 ← ovaj repo
  opengl-forward-renderer\     ← fork engina
```

Pokretanje iz foldera `PRS-Project` (PowerShell):

```powershell
cd ..\opengl-forward-renderer
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPRS_EROSION_DIR=../PRS-Project
cmake --build build --parallel
cd build\test-app
.\TestApp.exe
```

- `-DPRS_EROSION_DIR=../PRS-Project` znači da engine koristi **ovaj lokalni folder**, pa se svaka izmjena simulacije vidi odmah nakon `cmake --build`, bez commita. Bez ove opcije engine koristi submodul `external/prs-erosion`, tj. verziju PRS-Projecta fiksiranu u forku.
- Prvi build traje 3–4 minute (preuzimaju se i kompajliraju assimp, glfw, imgui...), a kasniji su brzi. Nakon prvog puta dovoljna su zadnja tri reda.
- `C:\msys64\ucrt64\bin` mora biti na PATH-u, kao za build ovog projekta.
- U aplikaciji: scena je na tipki **1**, kamera se pomjera s WASD i mišem, a **TAB** oslobađa miš za UI. Kontrole simulacije su u lijevom panelu, tab **"Terr"**.

Detalji (submodul, kontrole, šta se dešava u svakom frameu): [docs/erosion.md u forku](https://github.com/KerimSikalo/opengl-forward-renderer/blob/feature/erosion-scene/docs/erosion.md).

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
| [bench/](bench/) | Skripte za mjerenje i rezultati (`bench/results/*.csv`). |

Biblioteka `erosion` ne zavisi ni od kakve grafike, pa se kasnije može uključiti u [opengl-forward-renderer](https://github.com/bakir004/opengl-forward-renderer) kao zasebna scena.

## Odstupanja od Mei et al. 2007

| Šta | Mei et al. | Kod nas | Zašto |
|---|---|---|---|
| `Ks`, `Kd` | udio po koraku | brzina [1/s] × `dt` | rezultat ne zavisi od izbora `dt` |
| Kapacitet `C` | `Kc · sin α · \|v\|` | `× lmax(d)` (Šťava et al. 2008), `Kdmax = 0,5` | na ćelijama koje se isprazne u jednom koraku brzina je umjetno `l/dt` (vidi "Podešavanje parametara") |
| Transport sedimenta | semi-Lagrangian `s(x − v·Δt)` | **konzervativno kroz virtualne cijevi** (`--transport flux`, zadano); Mei-jeva metoda ostaje kao `--transport sl` | vidi ispod |

**Transport sedimenta.** Semi-Lagrangian advekcija rješava `∂s/∂t + v·∇s = 0`, dakle tretira vodu kao nestišljivu. U plitkoj vodi to ne vrijedi: tamo gdje se tok usporava i skuplja (dno padine, ulaz u jezero, `∇·v < 0`) sediment nestaje, a to su upravo mjesta gdje bi se trebao taložiti. Na fBm terenu 256² s kišom, nakon 5000 koraka:

| Metoda | Promjena ukupnog materijala | Promjena terena (min / max) |
|---|---|---|
| semi-Lagrangian (Mei) | **−19,8 %** | −70,3 / +31,6 |
| konzervativni kroz cijevi | **+1,2·10⁻⁷** | −29,4 / +21,2 |

### Podešavanje parametara: `Kdmax`

Formula za brzinu `u = ΔW / (l · d̄)` daje na ćelijama koje se isprazne u jednom koraku `|v| ≈ l / dt`. Mjerenje pokazuje da je to **artefakt, a ne fizika**: kad se `dt` prepolovi, udio mokrih ćelija s CFL > 0,5 skoro ne pada (za fizikalnu brzinu bi pao na pola), dakle brzina raste kao `1/dt`. Manji `dt` zato ne pomaže, a faktor dubine `lmax` (Šťava et al.) pomaže. fBm 256², kiša 0,02, 100 s simuliranog vremena; "jama" je ćelija niža za više od 1 od sva 4 susjeda:

| `Kdmax` | jame (dt = 0,02 / 0,01 / 0,005) | mokre ćelije s CFL > 0,5 | najdublja erozija |
|---|---|---|---|
| 0,1 | 191 / 258 / 259 | 26 % / 22 % / 16 % | −29 |
| **0,5** (zadano) | **31 / 39 / 28** | 50 % / 46 % / 39 % | −25 |
| 1,0 | 29 (dt = 0,02) | 60 % | −23 |

Uz `Kdmax = 0,5` erozija se skuplja u jasne kanale umjesto da se razlije po padini, a rezultat praktično ne zavisi od `dt`, što potvrđuje i skaliranje `Ks`/`Kd` s `dt`.

### Transport sedimenta

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

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
  - *Fast Hydraulic Erosion Simulation and Visualization on GPU* (Mei, Decaudin, Hu, PG 2007): **osnova algoritma** (virtual pipe model + erozija/depozicija). PDF treba ručno preuzeti, vidi [linkovi-implementacije.md](Izvori/linkovi-implementacije.md).
  - *Interactive Terrain Modeling Using Hydraulic Erosion* (Šťava et al., SCA 2008): proširenje na više slojeva materijala, odronjavanje obala i podjelu na blokove (tiling).
  - *A Second-Order Well-Balanced Positivity Preserving Central-Upwind Scheme for the Saint-Venant System* (Kurganov, Petrova 2007): numerički tačnija shema za shallow-water jednačine, kao poređenje s pipe modelom.
- **[Paralelizacije/](Izvori/Paralelizacije/)**: paralelne (GPU/CPU) implementacije shallow-water i erozije.
  - *Simulation and Visualization of the Saint-Venant System using GPUs* (Brodtkorb et al., 2010)
  - *Fast Hydraulic and Thermal Erosion on the GPU* (Jákó, CESCG 2011)
  - *3D Real-Time Hydraulic Erosion Simulation using Multi-Layered Heightmaps* (Nilles et al., VMV 2024)
  - Šira referenca (PDF ručno): *Efficient Shallow Water Simulations on GPUs* (Brodtkorb, Sætra, Altinakar 2012) i *TRITON: A Multi-GPU Open Source 2D Hydrodynamic Flood Model* (Morales-Hernández et al. 2021).
- **[linkovi-implementacije.md](Izvori/linkovi-implementacije.md)**: GitHub repozitoriji sa implementacijama, te linkovi za radove koje treba ručno preuzeti.

## Početni plan implementacije

C++20 i CMake. Teren je mreža `N × N` ćelija, a svaka ćelija čuva: visinu terena `b`, visinu vode `d`, sediment `s`, fluks prema 4 susjeda `f_L, f_R, f_T, f_B` i brzinu `(u, v)`. Podaci su u **SoA** rasporedu (posebni nizovi po veličini) zbog SIMD-a.

1. **Sekvencijalna verzija.** Jedan korak simulacije (po Mei et al.):
   1. dodavanje vode (kiša / izvor),
   2. računanje fluksa kroz virtualne cijevi,
   3. ažuriranje visine vode i polja brzine,
   4. erozija i depozicija (kapacitet nošenja sedimenta),
   5. transport sedimenta (semi-Lagrangian advekcija),
   6. isparavanje,
   7. (opcionalno) termalna erozija.

   Ulaz: generisani teren (Perlin/fBm šum) ili DEM. Izlaz: heightmap kao PGM/PNG.
   Metrike: vrijeme po koraku, broj ažuriranih ćelija u sekundi (cells/s), očuvanje mase vode i materijala.
2. **SIMD.** Kernele fluksa i erozije pišemo u tri verzije: skalarna, `#pragma omp simd` i ručni AVX2/FMA intrinsics.
3. **OpenMP.** Paralelizacija svakog koraka po redovima mreže (`parallel for`, poređenje `static` i `dynamic` rasporeda). Koristimo double-buffering da nema data race-a između čitanja susjeda i pisanja. Advekcija sedimenta je nepravilan pristup memoriji, pa nju posebno mjerimo. Probat ćemo i blokovsku (tiled) podjelu zbog keša.
4. **Mjerenja.** Veličine mreže ∈ {256², 512², 1024², 2048²}, niti ∈ {1, 2, 4, 8}, ubrzanje i efikasnost u odnosu na sekvencijalnu verziju, te strong i weak scaling.

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

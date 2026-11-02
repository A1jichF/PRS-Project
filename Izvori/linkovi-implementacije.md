# Linkovi za implementaciju

## GitHub repozitoriji (glavni)

Sva tri implementiraju **virtual pipe** shallow-water model iz rada Mei et al. 2007 (isti model koji planiramo mi), pa su direktno korisni kao referenca za strukturu koraka simulacije.

| # | Repozitorij | Jezik / platforma | Zašto je koristan |
|---|-------------|-------------------|-------------------|
| 1 | [simonmeister/hydraulic-terrain-modeler](https://github.com/simonmeister/hydraulic-terrain-modeler) | C++, CUDA, OpenGL | Najbliže našem projektu: C++ kod, simulacija vode i erozije paralelizovana u CUDA kernelima. |
| 2 | [bshishov/UnityTerrainErosionGPU](https://github.com/bshishov/UnityTerrainErosionGPU) | C#, Unity compute shaderi | Hidraulička + termalna erozija sa shallow-water jednačinama; vrlo čitljivo razdvojeni koraci (flux → voda → brzina → erozija/depozicija → transport → isparavanje). Dobar README s formulama. |
| 3 | [LanLou123/Webgl-Erosion](https://github.com/LanLou123/Webgl-Erosion) | TypeScript, WebGL | Interaktivna erozija u browseru (shaderi po koracima Mei et al.), proširena termalnom erozijom i drugim efektima; dobro za vizuelnu provjeru rezultata. |

## Dodatni repozitoriji (drugi pristupi, za poređenje)

- [SebLague/Hydraulic-Erosion](https://github.com/SebLague/Hydraulic-Erosion): erozija pomoću čestica (kapljica), a ne shallow-water mreže. Jednostavna sekvencijalna referenca.
- [dandrino/terrain-erosion-3-ways](https://github.com/dandrino/terrain-erosion-3-ways): Python, tri načina generisanja erodiranog terena (simulacija, ML, ...), s dobrim objašnjenjima.
- [Huw-man/Interactive-Erosion-Simulator-on-GPU](https://github.com/Huw-man/Interactive-Erosion-Simulator-on-GPU): C++ GPU simulator erozije.
- [TRITON (ORNL)](https://code.ornl.gov/hydro/triton): profesionalni multi-GPU/multi-CPU 2D shallow-water solver (C++, MPI, OpenMP, CUDA). Referenca za paralelizaciju, nije za eroziju.

## Ostalo

- Google Scholar pretraga: https://scholar.google.com/scholar?q=shallow+water+hydraulic+erosion+GPU
- Slajdovi uz Mei et al. 2007: https://slidetodoc.com/fast-hydraulic-erosion-simulation-and-visualization-on-gpu/

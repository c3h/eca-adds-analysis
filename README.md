# ECA ADDS Analysis

Supplementary code for the paper:

> **Correlating a measure of asynchrony of block-sequential updates with the dynamics of elementary cellular automata**
>
> D. Gimigliano, P.P. Balbi
> Universidade Presbiteriana Mackenzie
> ASCAT 2026

## Overview

This Wolfram Language code analyzes Elementary Cellular Automata (ECAs) under block-sequential update schedules (BSUs), computing:

- **Rotation-equivalent UDS representatives** for cyclic configurations
- **Asynchrony degree** ($\mathcal{T}$) for each BSU
- **Basins of attraction fields** for all (rule, BSU) pairs
- **Partition signatures** based on structural equivalence of basins
- **ADDS classes** (Asynchrony Degree based Dynamical Similarity)

## Requirements

- Wolfram Mathematica 12.0+ or Wolfram Engine
- Alternatively: [Wolfram Script](https://www.wolfram.com/wolframscript/)

## Usage

### In Mathematica

Open `ECA_ADDS_Analysis.wl` in Mathematica and evaluate all cells.

### From command line

```bash
wolframscript -file ECA_ADDS_Analysis.wl
```

### Parameters

Edit the parameters at the top of the file:

```mathematica
latticeSize = 5;  (* Configuration length L *)
radius = 1;       (* Neighbourhood radius *)
```

## Output

For `L=5`, the code identifies 10 ADDS classes among the 88 ECA equivalence classes, with partition signatures describing how BSUs with the same asynchrony degree produce structurally equivalent basins of attraction.

## Citation

If you use this code, please cite:

```bibtex
@inproceedings{gimigliano2026correlating,
  title={Correlating a measure of asynchrony of block-sequential updates
         with the dynamics of elementary cellular automata},
  author={Gimigliano, Daniel de Azevedo and Balbi, Pedro Paulo},
  booktitle={Asian Symposium on Cellular Automata Technology (ASCAT)},
  year={2026}
}
```

## License

MIT License - see [LICENSE](LICENSE) file.

## Acknowledgements

Based on CAMaT (Cellular Automata Mathematica Toolbox).
D. Gimigliano thanks Instituto Presbiteriano Mackenzie for a graduate grant.
P.P. Balbi is supported by CNPq under research grant PQ 303356/2022-7.

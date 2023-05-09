# ADinflation

Code computing the scalar spectral index, scalar-to-tensor ratio, number of e-folding, transfer functions for the two-field inflationary model based on the Afflec-Dine potential with non-minimal coupling to gravity. The physics underlying this model is discussed in

J. M. Cline, M. Puel, T. Toma, __Affleck-Dine inflation__, Phys. Rev. D 101, 043014 (2020), https://doi.org/10.1103/PhysRevD.101.043014

Please cite it (bibtex below) if you use this code.

## Structure
The main file is **ADinf.py**, which can be run in the usual way:
```
python ADinf.py
```
All the functions are defined in **myTools.py**. To run the Markov-Chain Monte-Carlo, use:
```
python MCMC.py
```
The plots are generated in **myPlots.py**.


### Bibtex citation
@article{Cline:2019fxx,
  author = "Cline, James M. and Puel, Matteo and Toma, Takashi",
  title = "{Affleck-Dine inflation}",
  eprint = "1909.12300",
  archivePrefix = "arXiv",
  primaryClass = "hep-ph",
  doi = "10.1103/PhysRevD.101.043014",
  journal = "Phys. Rev. D",
  volume = "101",
  number = "4",
  pages = "043014",
  year = "2020"
}


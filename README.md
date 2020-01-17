# exact-hubbard
Some exact solutions to the Hubbard model

This program calculates correlators in the Hubbard model.
The Hubbard model as used here is given by the hamiltonian[1]<br>
```
    H = -κ Σᵢⱼ (aᵢ^† aⱼ - bᵢ^† bⱼ) δ〈i,j〉 + U/2 Σᵢ (nᵢ - ñᵢ)²
```
See the paper for definitions of the symbols.

The single particle correlator is defined as
```
    Cᵢⱼ(τ) = 〈aᵢ(τ) aⱼ^†(0)〉 = Tr[aᵢ(τ) aⱼ^†(0) exp(-βH)] / Tr[exp(-βH)] 
```
where the trace is computed as 
```
    Tr[aᵢ^†(τ) aⱼ(0) exp(-βH)] 
        = Σₐₗₚₕₐ 〈alpha| aᵢ(τ) aⱼ^†(0) exp(-βH) |alpha〉
        = Σₐₗₚₕₐ Σ𝓰ₐₘₘₐ exp((τ-β) E_ₐₗₚₕₐ -τ H_𝓰ₐₘₘₐ) Aᵢᵃˡᵖʰᵃ ᵍᵃᵐᵐᵃ Aⱼᵍᵃᵐᵐᵃ ᵃˡᵖʰᵃ
```
with 
```
    Aᵢᵃˡᵖʰᵃ ᵍᵃᵐᵐᵃ = 〈alpha|aᵢ(0)|gamma〉.
```


## Conventions
There is some freedom in how creation and annihilation operators are defined.
This program uses the following for a state with a particle at site 0 and a hole at site 1
```
    |p   h〉 = a₀^† b₁^† |·  · 〉 
```  
and for a state with a particle and hole at site 0
```
    |ph  · 〉 = a₀^† b₀^† |·  · 〉 
```  


## Requirements
The main program needs
- C++17 compiler
- blaze
- LAPACK

Both need to be discoverable by CMake.

The analysis scripts require
- Python 3
- Numpy
- Matplotlib
 

## Usage
The main program can be compiled using CMake via
```shell script
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=RELEASE
cmake --build . -j <number-of-threads>
``` 
and run via the executable `exact_hubbard`.
This produces two files in the main directory:
- `spectrum.dat` contains the spectrum of the hamiltonian.
- `correlators.dat` contains the correlators.

There are rudimentary analysis / plot scripts written in Python in the `ana` directory.
They showcase how to read the data produced by `exact_hubbard`.  


## References
[1] Wynen, Berkowitz, Körber, Lähde, Luu
 *Avoiding Ergodicity Problems in Lattice Discretizations of the Hubbard Model*
 [DOI:10.1103/PhysRevB.100.075141](https://doi.org/10.1103/PhysRevB.100.075141)

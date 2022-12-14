# QuantumTransport

QuantumTransport is a package for calculating electron spin polarization and 
transmission of a molecular junction system. Calculations of electron 
transmission through a molecule (central region) connected to (semi-infinite)
contacts (left and right electrodes) are based on the non-self-consistent 
Green’s function method, an approximation to the non-equilibrium Green’s 
function (NEGF) method.

Spin polarization in the electron transmission of closed- and open shell
systems (such as radicals or diradicals) can be performed by interfacing with 
electronic structure packages. For closed-shell systems, it requires 
electronic structure calculations taking spin−orbit coupling (SOC) into 
account.  Greeen's functions of the left and right electrodes are approximated 
within the wide-band limit, using a constant local density of states, cLDOS.
Coupling terms of electrodes are neglected. 

## Requirements:
- Python3.x
- ase>=3.22.1
- cclib>=1.7.1
- matplotlib
- munch>=2.5.0
- numpy>=1.20.1

## Installation:
```
git clone https://github.com/nguyenka/QuantumTransport.git
cd QuantumTransport
pip --user install -e .
```

## Using qt script: 

The script for calculating electron spin polarization and transmission reads an output
file from Gaussian or Turbomole programs as an input file (filename):

```
$ qt filename [-ciss]
```

Results for spin polarization and transmission functions are written in 
transmission.txt.

Note that extra keywords must be added to a Gaussian job:
gfinput,iop(6/7=3),iop(5/33=3),iop(3/33=1).

For ```-ciss``` option with Turbomole the following keywords:
``` $intsdebug SAO ``` for scalar calculations.  
``` $soghf ```, ```$rx2c```, and ```$kramers``` for SOC calculations.  
File ```mos```, ```spinor.r```, and ```spinor.i``` are required in addition to the
output.

## Settings

The default settings for transport calculations. These settings can overwritten by
the user supplied in the `transport_params.json`. The keywords are given below.

|keyword | default | description|
|--------|---------|-----------|
|cLDOS | '0.036' | local density of states constant|
|Efermi | -5.0 | Fermi energy used for transport calculation in eV|
|eini | -2.0 | initial lower energy limit for the transmission function in eV|
|estep | 0.01| step size of the transmissions function in eV|
|eend | -7.0 | final energy limit for the transmission function in eV|

The settings are in json format:
{"cLDOS": 0.036, "Efermi": -5.0, "eini":-7.0, "eend": -2.0, "estep": 0.01}

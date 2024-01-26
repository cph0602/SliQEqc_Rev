# SliQEqc_Rev - A BDD-simulation-based Reversible Circuit Equivalence Checker

## Introduction
`SliQEqc_Rev` is a BDD-simulation-based reversible circuit equivalence checker implemented in C/C++ on top of [SliQSim]

## Build
To build the simulator, one needs to first configure `CUDD`:
```commandline
cd cudd
./configure --enable-dddmp --enable-obj --enable-shared --enable-static 
cd ..
```
Next, build the binary file `SliQSim`:
```commandline
make
```

## Execution
The circuit format being simulated is `OpenQASM` used by IBM's [Qiskit](https://github.com/Qiskit/qiskit), and the gate set supported in this simulator now contains Pauli-X (x), Controlled-NOT (cx), Toffoli (ccx and mcx). One can find some example benchmarks in [examples] folder. 

For example, to run the random benchmark 10_500.qasm in the [examples] folder, run:
```commandline
./SliQSim --sim_qasm examples/10_500.qasm
```

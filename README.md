# HOWTO

1. Define simulation prefix in file `run.sh` and parameters in `params.py`
2. send data to remote computer with `make send`
3.

# TODO

1. Dealiasing
2. Other analytical tests (Kloeden Platen p. 118)
3. Understand difference IF and ETD, Higher order

# SIMULATIONS

**Estelas**
0: JKW, alpha=0, dt=2 dx^2, works
1: Euler-Maruyama, alpha=0, dt=2 dx^2, works
2: JKW, alpha=0.1, dt=2 dx^2, seems to work?
3:
4: JKW, alpha=0.1, dt=2 dx^2, big parameters, works very well

**Tabua**
0: Euler-Maruyama, alpha=0, dt=2 dx^2, works
1: JKW, alpha=0, dt=2*dx^2
2: JKW, alpha=0, dt=0.2*dx^2

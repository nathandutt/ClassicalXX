# Classical Spin Chain Simulator

A script to simulate classical spin chains.

Originally designed to simulate the classical limit of the **XX spin chain model**, it has since been extended to **Landau-Lifshitz easy-axis dynamics** for studying quenches and coarsening — differing only in the effective magnetic field.

## Usage

Set simulation parameters in `simu.params`, then:

```bash
make          # compile
./ClassicalXX # run
```

To visualize results:

```bash
python3 plot.py  # press Space to advance in time
```

## TODO

- Add symplectic solver
- Rework evolution to reduce unnecessary copies of the chain

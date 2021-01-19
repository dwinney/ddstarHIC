# ddstarHIC
Investigations into X(3872) production through D-D* rescattering in Heavy Ion Collisions.

Most everything is done through the `data_file` class. This will take in an .dat file and convert it to .root format while calculating all relevant on-energies and missing mass and momentum in the lab and center-of-mass frames. 

A minimum running script can be executed with:
```bash
python data_file.py <path/to/.datfile>
```

Histograms and such are generated in the jupyter scripts.
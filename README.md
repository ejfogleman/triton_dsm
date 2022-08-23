# triton_dsm
## C-Simulator for "Triton" 2nd order DSADC (2nd order DSM w/ 2nd order mismatch-shaping DAC)
[Eric Fogleman, Jared Welz, Ian Galton, "An Audio ADC Delta–Sigma Modulator with 100-dB Peak SINAD and 102-dB DR Using a Second-Order Mismatch-Shaping DAC", IEEE JOURNAL OF SOLID-STATE CIRCUITS, VOL. 36, NO. 3, MARCH 2001.](https://ieeexplore.ieee.org/document/910472)  

The configuration file contains parameters corresponding to the published design.  (If I were to do this again today, I would make changes, but I left it as-is to reflect what was taped out and tested.)

Runs as a command-line program with the name of the input config file as its argument:

`> dsm t1_01`

Produces t1_01.dat containing output power spectral density data and t1_01.log with a copy of the run parameters.  The Matlab .m scripts or the Python library triton_dsm.py can be used to read and plot the program's output spectra.

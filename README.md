# wind-wave-directionality-paper
Codes and some data required to produce figures for a paper on wind wave directionality.

Manuscript currently in preparation for the *Journal of Physical Oceanography* (*JPO*):

"The directionality of wind waves from equilibrium to gravity-capillary scales"
Nathan J. M. Laxague, Z. Göksu Duvarci, Jan-Victor Björkqvist, Junzhe Liu, Lindsay Hogan, Alejandro Cifuentes-Lorenzen, and Christopher J. Zappa

My own MATLAB scripts/functions are licensed under [CC-BY](LICENSE).

## First Steps

Clone the repository:
```
git clone https://github.com/unh-cassll/wind-wave-directionality-paper.git
cd wind-wave-directionality-paper
```
Start MATLAB in your terminal:
```
matlab -nodisplay -nosplash
```

Call the data-grabbing routine (will run without input):
```
aa_step00_grab_data
```
Call the figure-production function.
This will run without arguments, but it takes two inputs:

1. list of figure numbers (array of integers, e.g., [1 3 5 7 9])
2. print figures to file? (Boolean, e.g. _true_ or _false_)

```
aa_step01_figure_generation_script([13 14],true) % e.g.
```
If you printed the figures to file (and you have access to a Linux CLI), you can convert the .svg files to .pdf by calling the bash script _svg2pdf_ from the terminal with the _figs/_ directory as an input:
```
./svg2pdf figs/
```

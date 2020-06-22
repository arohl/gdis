## GDIS, Copyright (C) 2000-2015 by Sean Fleming, Andrew Rohl

<andrew.rohl@curtin.edu.au>

GDIS comes with ABSOLUTELY NO WARRANTY.

This is free software. You are welcome to redistribute copies provided the conditions of the Version 2 GPL (GNU Public License) are met.

Although you are not required to do so, the authors would consider it a courtesy if you submit to them any changes you consider to be worthwhile. The goal would be to keep the development of GDIS more
or less centralized.

### Installation

To install GDIS from source you need a working C compiler, [gtk+2](http://www.gtk.org) and [gtkglext](https://projects.gnome.org/gtkglext/index.html). You then simply type:
```
./install
```
and follow the prompts. The script will check for the existence of gtk+2 and  gtkglext and then compile GDIS. If you don't specify an install directory, the GDIS executable will be in the `bin` directory.
 
### Notes

If you are compiling on a Mac, it is recommended that you use [MacPorts](https://www.macports.org) to install gtk and gtkglext. We have been unable to get it working with [Homebrew](https://brew.sh). If you know how to make it work with Homebrew, please get in touch!

Some cool features of GDIS are entirely due to the hard work that has gone into creating the following packages:

1. [CDD](http://www.inf.ethz.ch/personal/fukudak/cdd_home/) for the computation of halfspace intersections (morphology display).
	&copy;Komei Fukuda

2. [GPeriodic](http://www.frantz.fi/software/gperiodic.php) for pretty Periodic Table display and editing.
	&copy;1999 Kyle R. Burton 

3. [SgInfo](http://cci.lbl.gov/sginfo/) for Space Group generator lookup.
	&copy;1994-96 Ralf W. Grosse-Kunstleve

4. [Brute force symmetry analyzer](http://www.cobalt.chem.ucalgary.ca/ps/symmetry/).
	&copy;1996 S. Pachkovsky 

There are a few optional packages that enhance the GDIS experience.

1. For rendering (and subsequent image/movie viewing) you will need:
	- [POVRay](http://www.povray.org)
	- [ImageMagick](http://imagemagick.org)

2. Although GDIS supports output files from a large number of codes, input files can be created and run within GDIS for:
	- [GULP](http://nanochemistry.curtin.edu.au/gulp/)
	- [GAMESS](http://www.msg.chem.iastate.edu/GAMESS/)
	- [SIESTA](http://departments.icmab.es/leem/siesta/)
	- [Monty](http://www.vsc.science.ru.nl/deij/monty.html)
	- [VASP](http://www.vasp.at/)
	- [USPEX](http://www.uspex-team.org/en/uspex/overview)

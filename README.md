1 Content

Provide tools for getting usefull information from so-called NODE files (ie log files from ECMWF/Meteo France/ACCORD main binary MASTERODB)

2 Installation

Please download the project.

Tools can be used directly without any installation, but it may be better to install it.
Installation takes place by way of the Makefile. Its use needs 2 locations to be defined:
- binary files (ie shell scripts) are located by variable B (eg B=~/bin in this Makefile)
- side files (R scripts, HTML files, etc) are referenced by variable P (P=$(HOME)/proc in this Makefile)

Once adapted, just run make install.

3 Dependencies

The tools provided are based on shell scripts and R scripts. So, users need this software on their own machine.

Important thing, R scripts use an R package, mfnode, available on Github (see Github's project HaroldMTO/mfnode).

Please download the tarball preferentially. Then proceed to the desired installation by way of this command:
- R CMD INSTALL -l path/to/Rlibs mfnode_[version].tar.gz

You may then define one of R's environment variables (R_LIBS_USER or R_LIBS) pointing at its installation location.
- define this R variable: R_LIBS_USER=path/to/Rlibs
- or this other variable: R_LIBS=path/to/Rlibs
- run the shell script (that uses R software)

In case R should be loaded on user's machine, please look at the shell scripts and adapt (if necessary) the "module load" sequence for R.

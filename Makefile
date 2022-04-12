MAKEFLAGS += --no-print-directory

# ne pas mettre ~ pour P : il faut un chemin absolu
P = $(HOME)/proc/diags
B = ~/bin

.PHONY: build install diags

build:
	# rien

install:
	! git status --porcelain 2>/dev/null | grep -qvE "^\?\? "
	make diags
	make $B/epy_dump.py
	make $B/diag.sh
	make $B/spnorm.sh
	make $B/error.sh
	make $B/procmap.sh
	make $B/conf.sh
	if git status >/dev/null 2>&1; then \
		grep -q $(shell git log -1 --pretty=format:%h 2>/dev/null) $P/version || \
			git log -1 --oneline >> $P/version; \
	fi

diags:
	mkdir -p $P
	cp -pruv diag.R spnorm.R error.R procmap.R $P

$B/epy_dump.py: epy_dump.py
	cp -uv epy_dump.py $B

$B/diag.sh: diag.sh
	sed -re "s:diags=.+:diags=$P:" diag.sh > $B/diag.sh
	chmod a+x $B/diag.sh

$B/spnorm.sh: spnorm.sh
	sed -re "s:diags=.+:diags=$P:" spnorm.sh > $B/spnorm.sh
	chmod a+x $B/spnorm.sh

$B/error.sh: error.sh
	sed -re "s:diags=.+:diags=$P:" error.sh > $B/error.sh
	chmod a+x $B/error.sh

$B/procmap.sh: procmap.sh
	sed -re "s:diags=.+:diags=$P:" procmap.sh > $B/procmap.sh
	chmod a+x $B/procmap.sh

$B/conf.sh: conf.sh
	sed -re "s:diags=.+:diags=$P:" conf.sh > $B/conf.sh
	chmod a+x $B/conf.sh

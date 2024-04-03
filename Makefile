MAKEFLAGS += --no-print-directory

# ne pas mettre ~ pour P : il faut un chemin absolu
P = $(HOME)/proc/diagnode
B = ~/bin

.PHONY: build install diagnode

build:
	# rien

install:
	! git status --porcelain 2>/dev/null | grep -qvE "^\?\? "
	make diagnode
	make $B/error.sh
	make $B/conf.sh
	make $B/norms.sh
	make $B/setup.sh
	if git status >/dev/null 2>&1; then \
		grep -q $(shell git log -1 --pretty=format:%h 2>/dev/null) $P/version || \
			git log -1 --oneline >> $P/version; \
	fi

diagnode:
	mkdir -p $P
	cp -pruv norms.R error.R procmap.R setup.html norms.html $P

$B/error.sh: error.sh
	sed -re "s:diags=.+:diags=$P:" error.sh > $B/error.sh
	chmod a+x $B/error.sh

$B/conf.sh: conf.sh
	sed -re "s:diags=.+:diags=$P:" conf.sh > $B/conf.sh
	chmod a+x $B/conf.sh

$B/norms.sh: norms.sh
	sed -re "s:diags=.+:diags=$P:" norms.sh > $B/norms.sh
	chmod a+x $B/norms.sh

$B/setup.sh: setup.sh
	sed -re "s:diags=.+:diags=$P:" setup.sh > $B/setup.sh
	chmod a+x $B/setup.sh

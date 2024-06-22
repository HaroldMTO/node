MAKEFLAGS += --no-print-directory

# ne pas mettre ~ pour P : il faut un chemin absolu
P = $(HOME)/proc/node
B = ~/bin

.PHONY: build install node

build:
	# rien

install:
	! git status --porcelain 2>/dev/null | grep -qvE "^\?\? "
	make node
	make $B/norms.sh
	make $B/setup.sh
	if git status >/dev/null 2>&1; then \
		grep -q $(shell git log -1 --pretty=format:%h 2>/dev/null) $P/version || \
			git log -1 --oneline >> $P/version; \
	fi

node:
	mkdir -p $P
	cp -pruv gpnorms.R spnorms.R procmap.R runtime.R setup.html norms.html runtime.html $P

$B/norms.sh: norms.sh
	sed -re "s:node=.+:node=$P:" norms.sh > $B/norms.sh
	chmod a+x $B/norms.sh

$B/setup.sh: setup.sh
	sed -re "s:node=.+:node=$P:" setup.sh > $B/setup.sh
	chmod a+x $B/setup.sh

$B/runtime.sh: runtime.sh
	sed -re "s:node=.+:node=$P:" runtime.sh > $B/runtime.sh
	chmod a+x $B/runtime.sh

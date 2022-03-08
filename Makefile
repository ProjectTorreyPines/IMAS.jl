JULIA_PKG_DEVDIR ?= $(HOME)/.julia/dev
CURRENTDIR = $(shell pwd)

all:
	@echo 'IMAS makefile help'
	@echo ''
	@echo ' - make install  : install IMAS and all of its dependencies'
	@echo ' - make update   : update IMAS and all of its dependencies'
	@echo ''

install: install_IMAS
	julia -e '\
using Pkg;\
Pkg.activate();\
Pkg.develop(["IMAS", "IMASDD", "CoordinateConventions"]);\
Pkg.resolve();\
try Pkg.upgrade_manifest() catch end;\
'

install_IMAS: install_CoordinateConventions
	if [ ! -d "$(JULIA_PKG_DEVDIR)/IMAS" ]; then\
		julia -e 'using Pkg; Pkg.develop(url="git@github.com:ProjectTorreyPines/IMAS.jl.git");';\
	fi
	julia -e '\
using Pkg;\
Pkg.activate("$(JULIA_PKG_DEVDIR)/IMAS");\
Pkg.develop(["CoordinateConventions"]);\
Pkg.resolve();\
try Pkg.upgrade_manifest() catch end;\
'

install_IMASDD:
	if [ ! -d "$(JULIA_PKG_DEVDIR)/IMASDD" ]; then\
		julia -e 'using Pkg; Pkg.develop(url="git@github.com:ProjectTorreyPines/IMASDD.jl.git");';\
	fi
	julia -e '\
using Pkg;\
Pkg.activate("$(JULIA_PKG_DEVDIR)/IMASDD");\
Pkg.resolve();\
try Pkg.upgrade_manifest() catch end;\
'

install_CoordinateConventions:
	if [ ! -d "$(JULIA_PKG_DEVDIR)/CoordinateConventions" ]; then\
		julia -e 'using Pkg; Pkg.develop(url="git@github.com:ProjectTorreyPines/CoordinateConventions.jl.git");';\
	fi
	julia -e '\
using Pkg;\
Pkg.activate("$(JULIA_PKG_DEVDIR)/CoordinateConventions");\
Pkg.resolve();\
try Pkg.upgrade_manifest() catch end;\
'

update: update_IMAS update_IMASDD update_CoordinateConventions
	make install

update_IMAS:
	cd $(JULIA_PKG_DEVDIR)/IMAS; git fetch; git pull

update_IMASDD:
	cd $(JULIA_PKG_DEVDIR)/IMASDD; git fetch; git pull

update_CoordinateConventions:
	cd $(JULIA_PKG_DEVDIR)/CoordinateConventions; git fetch; git pull

.PHONY:

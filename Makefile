# ---- toolchain ------------------------
CC       ?= cc
CFLAGS   ?= -O3 -std=c11 -Wall -Wextra
CPPFLAGS ?= -Ibin/seqtk
LDLIBS   ?= -lz

# CORES used by snakemake for LG_db builds. Override: `make CORES=8`.
CORES    ?= $(shell nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2)


all: BCnS_ALG UnicellLG UnicellLGOrthofinder BCnSSimakov2022 seqtk fainfo

.PHONY: seqtk fainfo

# ---- fainfo -----------------------------------------------------------------
# Build ./bin/fainfo from ./source/fainfo.c using kseq.h in bin/seqtk/
fainfo: bin/fainfo

KSEQ_URL = https://raw.githubusercontent.com/lh3/seqtk/master/kseq.h

bin/fainfo: source/fainfo.c bin/seqtk/seqtk
	@mkdir -p bin
	@curl -L $(KSEQ_URL) -o source/kseq.h
	@$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(LDLIBS)
	@echo "Built: $@"

seqtk: bin/seqtk/seqtk

bin/seqtk/seqtk:
	@mkdir -p bin
	@if [ ! -d bin/seqtk ]; then \
	  git clone https://github.com/lh3/seqtk.git bin/seqtk; \
	else \
	  echo "bin/seqtk exists, assuming repo already cloned."; \
	fi
	@echo "Building seqtk..."
	@$(MAKE) -C bin/seqtk
	@echo "Built: $@"

BCnS_ALG: LG_db/BCnS_LGs.tar.gz
	cd LG_db; \
	tar --skip-old-files -xzvf BCnS_LGs.tar.gz; \
	cd BCnS_LGs; \
	snakemake --cores $(CORES)

BCnSSimakov2022: LG_db/BCnSSimakov2022.tar.gz
	cd LG_db; \
	tar --skip-old-files -xzvf BCnSSimakov2022.tar.gz; \
	cd BCnSSimakov2022; \
	snakemake --cores $(CORES)

UnicellLG: LG_db/UnicellMetazoanLgs.tar.gz
	cd LG_db; \
	tar --skip-old-files -xzvf UnicellMetazoanLgs.tar.gz; \
	cd UnicellMetazoanLgs; \
	snakemake --cores $(CORES)

UnicellLGOrthofinder: LG_db/UnicellMetazoanLgsOrthofinder.tar.gz
	cd LG_db; \
	tar --skip-old-files -xzvf UnicellMetazoanLgsOrthofinder.tar.gz; \
	cd UnicellMetazoanLgsOrthofinder; \
	snakemake --cores $(CORES)

CLG_v1.0: LG_db/CLG_v1.0.tar.gz
	cd LG_db; \
	tar --skip-old-files -xzvf CLG_v1.0.tar.gz; \
	cd CLG_v1.0; \
	snakemake --cores $(CORES)

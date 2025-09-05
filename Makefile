all: BCnS_ALG UnicellLG UnicellLGOrthofinder seqtk

.PHONY: seqtk

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
	snakemake

UnicellLG: LG_db/UnicellMetazoanLgs.tar.gz
	cd LG_db; \
	tar --skip-old-files -xzvf UnicellMetazoanLgs.tar.gz; \
	cd UnicellMetazoanLgs; \
	snakemake

UnicellLGOrthofinder: LG_db/UnicellMetazoanLgsOrthofinder.tar.gz
	cd LG_db; \
	tar --skip-old-files -xzvf UnicellMetazoanLgsOrthofinder.tar.gz; \
	cd UnicellMetazoanLgsOrthofinder; \
	snakemake

CLG_v1.0: LG_db/CLG_v1.0.tar.gz
	cd LG_db; \
	tar --skip-old-files -xzvf CLG_v1.0.tar.gz; \
	cd CLG_v1.0; \
	snakemake

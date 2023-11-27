CORES=2
all: BCnS_ALG UnicellLG UnicellLGOrthofinder BCnSSimakov2022

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

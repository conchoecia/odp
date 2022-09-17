all: BCnS_ALG UnicellLG UnicellLGOrthofinder

BCnS_ALG: LG_db/BCnS_LGs.tar.gz
	cd LG_db; \
	tar -xzvf BCnS_LGs.tar.gz; \
	cd BCnS_LGs; \
	snakemake

UnicellLG: LG_db/UnicellMetazoanLgs.tar.gz
	cd LG_db; \
	tar -xzvf UnicellMetazoanLgs.tar.gz; \
	cd UnicellMetazoanLgs; \
	snakemake

UnicellLGOrthofinder: LG_db/UnicellMetazoanLgsOrthofinder.tar.gz
	cd LG_db; \
	tar -xzvf UnicellMetazoanLgsOrthofinder.tar.gz; \
	cd UnicellMetazoanLgsOrthofinder; \
	snakemake

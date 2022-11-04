all: BCnS_ALG UnicellLG UnicellLGOrthofinder

BCnS_ALG: LG_db/BCnS_LGs.tar.gz
	cd LG_db; \
	tar -xzvf BCnS_LGs.tar.gz; \
	cd BCnS_LGs; \
	snakemake --cores 1

UnicellLG: LG_db/UnicellMetazoanLgs.tar.gz
	cd LG_db; \
	tar -xzvf UnicellMetazoanLgs.tar.gz; \
	cd UnicellMetazoanLgs; \
	snakemake --cores 1 

UnicellLGOrthofinder: LG_db/UnicellMetazoanLgsOrthofinder.tar.gz
	cd LG_db; \
	tar -xzvf UnicellMetazoanLgsOrthofinder.tar.gz; \
	cd UnicellMetazoanLgsOrthofinder; \
	snakemake --cores 1

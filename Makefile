all: BCnS_ALG

BCnS_ALG: LG_db/BCnS_LGs.tar.gz
	cd LG_db; \
	tar -xzvf BCnS_LGs.tar.gz; \
	cd BCnS_LGs; \
	snakemake

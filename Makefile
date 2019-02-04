
# $^ == all the dependencies
# $@ == the target
# within make $ means variable
# $^, $@, etc special variables
# $(...) your declared variables

DATADIR ?= ~/cancerData

ALLSRCS := $(shell cd $(DATADIR); ls *.csv)
ALLRDS := $(ALLSRCS:csv=rds)

default: $(ALLRDS)

.PHONY: test

test:
	@echo $(ALLRDS)

%.rds: data_preprocessing.R $(DATADIR)/%.csv
	Rscript $^ $@

cleanrds:
	rm $(ALLRDS)
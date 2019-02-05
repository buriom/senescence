
# $^ == all the dependencies
# $@ == the target
# within make $ means variable
# $^, $@, etc special variables
# $(...) your declared variables

DATADIR ?= ~/cancerData
#MODELDIR?= ~/senescenceModel
ALLSRCS := $(shell cd $(DATADIR); ls *.csv)
ALLRDS := $(ALLSRCS:csv=rds)

default: $(ALLRDS)

.PHONY: test

test:
	@echo $(ALLRDS)

%.rds: data_preprocessing.R $(DATADIR)/%.csv
	Rscript $^ $@

#%.rds fitting.R $(MODELDIR)/model%.Rds
cleanrds:
	rm $(ALLRDS)
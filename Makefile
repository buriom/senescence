default: alldatasrcs

test:
	@echo $(FIGJPG)
# $^ == all the dependencies
# $@ == the target
# within make $ means variable
# $^, $@, etc special variables
# $(...) your declared variables

#Data processing
DATADIR ?= cancerData
PROCESSEDIR ?= preProcessed
FITSDIR ?= fits
FIGDIR ?= figures

ALLSRCS := $(shell cd $(DATADIR); ls *.csv)
ALLRDS := $(addprefix $(PROCESSEDIR)/,$(ALLSRCS:csv=rds))
FITSRDS := $(addprefix $(FITSDIR)/fit_simple_,$(ALLSRCS:csv=rds))
FIGJPG := $(addprefix $(FIGDIR)/,$(ALLSRCS:csv=jpg))

alldatasrcs: $(FIGJPG)
	
$(PROCESSEDIR):
	mkdir $@
 
$(PROCESSEDIR)/%.rds: data_preprocessing.R $(DATADIR)/%.csv | $(PROCESSEDIR)  
	Rscript $^ $@


#model specification
model_%.rds: senescence_model_%.R 
	Rscript $^ $@

#model fitting
$(FITSDIR):
	mkdir $@
	
$(FITSDIR)/fit_simple_%.rds: fitting.R model_simple.rds $(PROCESSEDIR)/%.rds | $(FITSDIR)
	Rscript $^ $@

$(FIGDIR):
	mkdir $@
	
$(FIGDIR)/%.jpg: plotting.R $(FITSDIR)/fit_simple_%.rds model_simple.rds $(PROCESSEDIR)/%.rds | $(FIGDIR)
	Rscript $^ $@

clean:
	rm $(ALLRDS)
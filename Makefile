default: alldatasrcs

test:
	@echo $(ALLDWNLDS)
# $^ == all the dependencies
# $@ == the target
# within make $ means variable
# $^, $@, etc special variables
# $(...) your declared variables
# $^ == all the dependencies
# $@ == the target


-include local.mk # optional way to define DATADIR, PROCESSEDIR, etc

# project file structure
DWNLDSDIR ?= cancerDownloads
DATADIR ?= cancerData
PROCESSEDIR ?= preProcessed
FITSDIR ?= fits
FIGDIR ?= figures
FIGFORMAT ?= jpg

$(DWNLDSDIR) $(DATADIR) $(PROCESSEDIR) $(FITSDIR) $(FIGDIR):
	mkdir $@


# targets for obtaining raw data






ALLSRCS := $(shell cd $(DATADIR); ls *.csv)
ALLDWNLDS := $(shell cd $(DWNLDSDIR); ls *.csv)
ALLRDS := $(addprefix $(PROCESSEDIR)/,$(ALLSRCS:csv=rds))
FITSRDS := $(addprefix $(FITSDIR)/fit_simple_,$(ALLSRCS:csv=rds))
FIGJPG := $(addprefix $(FIGDIR)/,$(ALLSRCS:csv=jpg))

alldatasrcs: $(FIGJPG)

#extractData:
#	for i in ALLDWNLDS; do 	cut -d "," -f 3-4  i | head -n 23 |tail -n 3+ > $(cut -d "," -f 1 i | head -n 1 .csv 

	
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
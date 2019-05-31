comma := ,
empty :=
  space := $(empty) $(empty)

default: test

test:
	@echo '$(FIGNAMES)'
# $^ == the dependencies before '|'
# $| == the dependencies after '|'
# $@ == the target
# within make $ means variable
# $^, $@, etc special variables
# $(...) your declared variables
# $^ == all the dependencies
# $@ == the target

# project file structure

PROCESSEDIR ?= preProcessed
#FIGDIR ?= figures
#FIGFORMAT ?= jpg
FIG1DIR ?= figures/model1Figs
FIG2DIR ?= figures/model2Figs
FIG3DIR ?= figures/model3Figs
#FITSDIR ?= fits
FITS1DIR ?= fits/model1fits
FITS2DIR ?= fits/model2fits
FITS3DIR ?= fits/model3fits
$(PROCESSEDIR) $(FIG1DIR) $(FIG2DIR) $(FIG3DIR) $(FITS1DIR) $(FITS2DIR) $(FITS3DIR):
	mkdir $@
  
# this means that make will look for these targets / dependencies in
# these directories

vpath %_data.rds $(PROCESSEDIR)

#plotting figures
ALLSRCS := $(shell cd $(PROCESSEDIR); ls *_data.rds)
#ALLSRCS := ProstateCancer_proc.csv TestisCancer_proc.csv BreastCancer_proc.csv CervixUteriCancer_proc.csv CervixUteriCancer_proc.csv CorpusandUterusNOSCancer_proc.csv FemalegenitalsystemCancer_proc.csv OvaryCancer_proc.csv VaginaCancer_proc.csv VulvaCancer_proc.csv

FIG1JPG := $(addprefix $(FIG1DIR)/,$(ALLSRCS:_data.rds=.jpg))
FIG2JPG := $(addprefix $(FIG2DIR)/,$(ALLSRCS:_data.rds=.jpg))
FIG3JPG := $(addprefix $(FIG3DIR)/,$(ALLSRCS:_data.rds=.jpg))
FIGNAMES := $(FIG1JPG) $(FIG2JPG)# $(FIG3JPG)


FITS1RDS := $(addprefix $(FITS1DIR)/,$(ALLSRCS:_data.rds=.rds))
FITS2RDS := $(addprefix $(FITS2DIR)/,$(ALLSRCS:_data.rds=.rds))
FITS3RDS := $(addprefix $(FITS3DIR)/,$(ALLSRCS:_data.rds=.rds))
FITSNAMES := $(FITS1RDS) $(FITS2RDS) $(FITS3RDS)


allfigs: $(FIG3JPG)
#	@echo $(FIGNAMES)

$(FIG1DIR)/%.jpg: model1.R $(PROCESSEDIR)/%_data.rds | $(FIG1DIR)
	-Rscript $^ $@
  
$(FIG2DIR)/%.jpg: model2.R $(PROCESSEDIR)/%_data.rds | $(FIG2DIR)
	-Rscript $^ $@
  
$(FIG3DIR)/%.jpg: model3Annual.R $(PROCESSEDIR)/%_data.rds | $(FIG3DIR)
	-Rscript $^ $@
  
  
  #fits
$(FITS1DIR)/%.rds: model1.R $(PROCESSEDIR)/%_data.rds | $(FITS1DIR)
	-Rscript $^ $@
  
$(FITS2DIR)/%.rds: model2.R $(PROCESSEDIR)/%_data.rds | $(FITS2DIR)
	-Rscript $^ $@
  
$(FITS3DIR)/%.rds: model3.R $(PROCESSEDIR)/%_data.rds | $(FITS3DIR)
	-Rscript $^ $@


cleanrds:
	rm -rf $(PROCESSEDIR)

cleanfig1:
	rm -rf $(FIG1DIR)

cleanfig2:
	rm -rf $(FIG2DIR)

cleanfig3:
	rm -rf $(FIG3DIR)
# ALLSRCS := $(shell cd $(DATADIR); ls *.csv)
# ALLDWNLDS := $(shell cd $(DWNLDSDIR); ls *.csv)
# ALLRDS := $(addprefix $(PROCESSEDIR)/,$(ALLSRCS:csv=rds))
# FITSRDS := $(addprefix $(FITSDIR)/fit_simple_,$(ALLSRCS:csv=rds))
# FIGJPG := $(addprefix $(FIGDIR)/,$(ALLSRCS:csv=jpg))
# 
# alldatasrcs: $(FIGJPG)
# 
# #extractData:
# #	for i in ALLDWNLDS; do 	cut -d "," -f 3-4  i | head -n 23 |tail -n 3+ > $(cut -d "," -f 1 i | head -n 1 .csv 
# 
# 
# $(PROCESSEDIR):
# 	mkdir $@
# 
# $(PROCESSEDIR)/%.rds: dataPreprocessing.R $(DATADIR)/%.csv | $(PROCESSEDIR)  
# 	Rscript $^ $@
# 
# 
# #model specification
# model_%.rds: senescence_model_%.R 
# 	Rscript $^ $@
# 
# #model fitting
# $(FITSDIR):
# 	mkdir $@
# 
# $(FITSDIR)/fit_simple_%.rds: fitting.R model_simple.rds $(PROCESSEDIR)/%.rds | $(FITSDIR)
# 	Rscript $^ $@
# 
# $(FIGDIR):
# 	mkdir $@
# 
# $(FIGDIR)/%.jpg: plotting.R $(FITSDIR)/fit_simple_%.rds model_simple.rds $(PROCESSEDIR)/%.rds | $(FIGDIR)
# 	Rscript $^ $@
# 
# clean:
# 	rm $(ALLRDS)
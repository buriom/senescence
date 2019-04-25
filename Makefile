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


-include local.mk # optional way to define DATADIR, PROCESSEDIR, etc

# project file structure
DWNLDSDIR ?= cancerDownloads
DATADIR ?= cancerData
PROCESSEDIR ?= preProcessed
#FITSDIR ?= fits
#FIGDIR ?= figures
FIGFORMAT ?= jpg
FIG1DIR ?= figures/model1Figs
FIG2DIR ?= figures/model2Figs
FIG3DIR ?= figures/model3Figs

$(DWNLDSDIR) $(DATADIR) $(PROCESSEDIR) $(FIG1DIR) $(FIG2DIR) $(FIG31DIR):
	mkdir $@

# targets for obtaining raw data

DATAURL := https://seer.cancer.gov/explorer/download_data.php
DATAQPARS := data_type=1&graph_type=3&compareBy=race&chk_sex_1=1&chk_race_1=1&hdn_data_type=1&advopt_precision=1&showDataFor=sex_1
#DATAIDS ?= 34 50 76 500
DATAIDS ?= 34 50 76 500 55 57 402 20 21 31 58 16 17 75 56 8 38 9 83 110 72 46 90 91 92 93 95 96 97 100 418 5 35 47 53 111 409 89 86 3 12 61 40 66 7 19 51 18 67 80 6 71 62 63

allrawdata: $(addsuffix _raw.csv,$(DATAIDS))

# this means that make will look for these targets / dependencies in
# these directories
vpath %_raw.csv $(DWNLDSDIR)
vpath %_proc.csv $(DATADIR)
vpath %_data.rds $(PROCESSEDIR)

%_raw.csv: | $(DWNLDSDIR)
	cd $|; wget -O $@ $(DATAURL)?site=$*&$(DATAQPARS)

## COMPLICATED PAIRING OF NAME_proc.csv TO ID_raw.csv
AVAILDLS := $(shell cd $(DWNLDSDIR); ls *_raw.csv)
findname = $(shell head -n 1 $(DWNLDSDIR)/$(id) | sed -e 's/^[[:blank:]]*//' -e 's/[[:blank:]]*$$//' -e 's/^"//' -e 's/"$$//' -e 's/[^[:alpha:]]//g')
PROCNAMES := $(foreach id,$(AVAILDLS),$(findname)_proc.csv)
PAIRS := $(join $(PROCNAMES),$(patsubst %,^%,$(AVAILDLS)))

#I Substituted tail -r with tac
define NAME_template
$(subst ^,: ,$(pr)) | $(DATADIR) # this specifies NAME_proc.csv: ID_raw.csv | DATADIR
	tail -n +4 $$^ | tac | tail -n +15 | tac | csvtool -t , col 1-7 - > $$|/$$@

endef

$(foreach pr,$(PAIRS),$(eval $(NAME_template))) #$(subst ^,$(comma),$(pr)))

allprocdata: $(PROCNAMES) 

RDSNAMES := $(PROCNAMES:_proc.csv=_data.rds)

allrdsreads: $(RDSNAMES)


%_data.rds: dataPreprocessing.R %_proc.csv | $(PROCESSEDIR)
	Rscript $^ $|/$@

#plotting figures
ALLSRCS := $(shell cd $(DATADIR); ls *_proc.csv)
FIG1JPG := $(addprefix $(FIG1DIR)/,$(ALLSRCS:_proc.csv=.jpg))
FIG2JPG := $(addprefix $(FIG2DIR)/,$(ALLSRCS:_proc.csv=.jpg))
FIG3JPG := $(addprefix $(FIG3DIR)/,$(ALLSRCS:_proc.csv=.jpg))

FIGNAMES := $(FIG1JPG) $(FIG2JPG) $(FIG3JPG)

allfigs: $(FIGNAMES)
#	@echo $(FIG1DIR)

$(FIG1DIR)/%.jpg: model1.R $(PROCESSEDIR)/%_data.rds | $(FIG1DIR)
	Rscript $^ $@

$(FIG2DIR)/%.jpg: model2.R $(PROCESSEDIR)/%_data.rds | $(FIG2DIR)
	Rscript $^ $@

$(FIG3DIR)/%.jpg: model3.R $(PROCESSEDIR)/%_data.rds | $(FIG3DIR)
	Rscript $^ $@
	
#cleaning
cleandl:
	rm -rf $(DWNLDSDIR)

cleanproc:
	rm -rf $(DATADIR)

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
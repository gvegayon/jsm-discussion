all:
	Rscript --vanilla R/prototype.R > R/prototype.Rout 2>&1

prepare:
	rm -rf package/*; \
	# Check if a directory exists
	if [ ! -d "package" ]; then mkdir package; fi; \
	cp discussion.tex bibliography.bib notation.tex package; \
	cp figures/max-out-stats.pdf package/.; \
	cp figures/conditional-prob-ttriad.pdf package/.; \
	cp figures/vif-n=5.pdf package/.; \
	# Copying R scripts (check if folder exists)
	if [ ! -d "package/R" ]; then mkdir package/R; fi; \
	# Copying R scripts into a single zip file
	zip -9 package/R.zip R/*R; \
	cp discussion.pdf package/.
	

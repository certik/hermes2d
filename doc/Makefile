# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         =

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d _build/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

.PHONY: help clean html web pickle htmlhelp latex changes linkcheck

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html      to make standalone HTML files"
	@echo "  pickle    to make pickle files"
	@echo "  json      to make JSON files"
	@echo "  htmlhelp  to make HTML files and a HTML help project"
	@echo "  latex     to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  changes   to make an overview over all changed/added/deprecated items"
	@echo "  linkcheck to check all external links for integrity"

clean:
	-rm -rf _build/*

html:
	mkdir -p _static
	mkdir -p _build/html _build/doctrees
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) _build/html
	@echo
	@echo "Build finished. The HTML pages are in _build/html."

pickle:
	mkdir -p _build/pickle _build/doctrees
	$(SPHINXBUILD) -b pickle $(ALLSPHINXOPTS) _build/pickle
	@echo
	@echo "Build finished; now you can process the pickle files."

web: pickle

json:
	mkdir -p _build/json _build/doctrees
	$(SPHINXBUILD) -b json $(ALLSPHINXOPTS) _build/json
	@echo
	@echo "Build finished; now you can process the JSON files."

htmlhelp:
	mkdir -p _build/htmlhelp _build/doctrees
	$(SPHINXBUILD) -b htmlhelp $(ALLSPHINXOPTS) _build/htmlhelp
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in _build/htmlhelp."

latex:
	mkdir -p _build/latex _build/doctrees
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) _build/latex
	@echo
	@echo "Build finished; the LaTeX files are in _build/latex."
	@echo "Run \`make all-pdf' or \`make all-ps' in that directory to" \
	      "run these through (pdf)latex."

changes:
	mkdir -p _build/changes _build/doctrees
	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) _build/changes
	@echo
	@echo "The overview file is in _build/changes."

linkcheck:
	mkdir -p _build/linkcheck _build/doctrees
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) _build/linkcheck
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in _build/linkcheck/output.txt."

# Below is the build script for .tex docs:

MAIN = manual.tex

IMGDIR = img


#### And ignore the rest! :-) ##################################################

PDFNAME = $(MAIN:.tex=.pdf)
IDXNAME = $(MAIN:.tex=.idx)

PSIMG  = $(wildcard $(IMGDIR)/*.ps)
EPSIMG = $(wildcard $(IMGDIR)/*.eps)
SVGIMG = $(wildcard $(IMGDIR)/*.svg)

PSPDF  = $(PSIMG:.ps=.pdf)
EPSPDF = $(EPSIMG:.eps=.pdf)
SVGPDF = $(SVGIMG:.svg=.pdf)


# always recompile the main PDF
.PHONY: $(PDFNAME) clean-pdf showpdf

$(PDFNAME): $(PSPDF) $(EPSPDF) $(SVGPDF)
	@-if which remark 2> /dev/null > /dev/null; then \
		pdflatex --interaction nonstopmode $(MAIN) | remark colors; \
	else \
		pdflatex --interaction nonstopmode $(MAIN); \
	fi
	@if [ -f $(IDXNAME) ]; then echo; echo; makeindex $(basename $(MAIN)); fi


# turn off GhostScript automatic page rotation
GSNOROTATE = GS_OPTIONS="-dAutoRotatePages=/None"

$(IMGDIR)/%.pdf: $(IMGDIR)/%.ps
	$(GSNOROTATE) epstopdf $<

$(IMGDIR)/%.pdf: $(IMGDIR)/%.eps
	$(GSNOROTATE) epstopdf $<

$(IMGDIR)/%.pdf: $(IMGDIR)/%.svg
	TMPFILE=`mktemp`; inkscape --file=$< --export-eps=$$TMPFILE --export-text-to-path; \
	$(GSNOROTATE) epstopdf $$TMPFILE --outfile=$@; rm $$TMPFILE


# remove all generated files (except the main PDF)
clean-pdf:
	@rm -f *.log $(PSPDF) $(EPSPDF) $(SVGPDF)

# show the names of all generated PDF files
showpdf:
	@echo $(PSPDF) $(EPSPDF) $(SVGPDF)


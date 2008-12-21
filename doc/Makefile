################################################################################
##
##  This is a generic pdfLaTeX makefile. Place it in the document directory and
##  change the main .tex file name and image directory name below. Then just
##  type 'make' to have your document built.
##
##  This file automatically converts all .ps, .eps and .svg images in your image
##  directory to their .pdf versions to be included in your pdfLaTeX document.
##
##  Example:
##
##  1. Create/modify a figure and save it as img/figure1.svg
##  2. In your document, write \includegraphics{img/figure1}
##  3. Type 'make'. The figure is automatically converted to img/figure1.pdf
##     and your document is compiled with pdflatex.
##
##  The same applies for other figure formats (.ps, .eps).
##
##  .jpg, .png and .pdf figures can be placed in the image directory too, but
##  since they need no processing, they are left alone.
##
##  The .svg conversion assumes you have inkscape installed.
##
##  To remove all automatically generated files, type 'make clean'.
##

#### Just change the following: ################################################

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
.PHONY: $(PDFNAME) clean showpdf

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
clean:
	@rm -f *.log $(PSPDF) $(EPSPDF) $(SVGPDF)

# show the names of all generated PDF files
showpdf:
	@echo $(PSPDF) $(EPSPDF) $(SVGPDF)


latexfile = template
TEX = latex

# support subdirectories
VPATH = Figs

# reruns latex if needed.  to get rid of this capability, delete the
# three lines after the rule.  Delete .bbl dependency if not using
# BibTeX references.
# idea from http://ctan.unsw.edu.au/help/uk-tex-faq/Makefile
$(latexfile).dvi : $(latexfile).tex
	$(TEX) $(latexfile)

# keep .eps files in the same directory as the .fig files

$(latexfile).pdf : $(latexfile).dvi
	dvipdf $(latexfile).dvi $(latexfile).pdf

pdf : $(latexfile).pdf



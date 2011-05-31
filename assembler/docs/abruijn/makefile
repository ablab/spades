default: abruijn

fig1.pdf: fig1.dot
	dot -Tpdf -o $@ $<

fig2.pdf: fig2.dot
	dot -Tpdf -o $@ $<

fig3.pdf: fig3.dot
	dot -Tpdf -o $@ $<

fig4.pdf: fig4.dot
	dot -Tpdf -o $@ $<

TAAACGAAAC.pdf: TAAACGAAAC.dot
	dot -Tpdf -o $@ $<

ATGCATTGCACTGCA.pdf: ATGCATTGCACTGCA.dot
	dot -Tpdf -o $@ $<

ATGCATTGCACTGCA_debruijn.pdf: ATGCATTGCACTGCA_debruijn.dot
	dot -Tpdf -o $@ $<

ATGCATTGCACTGCA_earmark.pdf: ATGCATTGCACTGCA_earmark.dot
	dot -Tpdf -o $@ $<

clean:
	rm *.log *.aux *.out *.toc *.tdo abruijn.pdf fig1.pdf fig2.pdf fig3.pdf fig4.pdf TAAACGAAAC.pdf mccmb_poster.pdf ATGCATTGCACTGCA.pdf ATGCATTGCACTGCA_debruijn.pdf ATGCATTGCACTGCA_earmark.pdf ATGCATTGCACTGCA_earmark2.pdf

pics: fig1.pdf fig2.pdf fig3.pdf fig4.pdf TAAACGAAAC.pdf ATGCATTGCACTGCA.pdf ATGCATTGCACTGCA_debruijn.pdf ATGCATTGCACTGCA_earmark.pdf 

abruijn: pics 
	pdflatex abruijn
	bibtex abruijn

poster: pics 
	pdflatex mccmb_poster

all: abruijn poster

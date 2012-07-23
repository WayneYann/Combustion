
default:
	latex SMCnotes
	latex SMCnotes
	dvips -t letter -o SMCnotes.ps SMCnotes.dvi
	ps2pdf -dPDFSETTINGS=/prepress SMCnotes.ps
clean:
	rm -rf *.aux *.dvi *.log *.toc *.lof *.lot *.blg *.bbl *~
realclean:
	rm -rf *.aux *.dvi *.log *.toc *.lof *.lot *.blg *.bbl *~ SMCnotes.ps 

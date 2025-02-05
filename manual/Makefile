.PHONY: clean html

html: ebmb.pdf
	sphinx-build -M html . .

ebmb.pdf: ebmb.tex ../logo/ebmb.pdf
	pdflatex --interaction=batchmode $<
	pdflatex --interaction=batchmode $<

../logo/ebmb.pdf: ../logo/ebmb.tex
	$(MAKE) -C ../logo/ ebmb.pdf

clean:
	rm -rf doctrees html

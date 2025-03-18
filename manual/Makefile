.PHONY: html examples clean

html: ebmb.pdf examples
	sphinx-build -M html . .

ebmb.pdf: ebmb.tex ../logo/ebmb.pdf
	pdflatex --interaction=batchmode $<
	pdflatex --interaction=batchmode $<

../logo/ebmb.pdf: ../logo/ebmb.tex
	$(MAKE) -C ../logo/ ebmb.pdf

examples:
	$(MAKE) -C ../examples/

clean:
	rm -rf doctrees html

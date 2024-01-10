.PHONY: all clean

all: ebmb.pdf ebmb.png ebmb_square.png ebmb_banner.png

%.pdf: %.tex
	pdflatex --interaction=batchmode $*

ebmb.png: ebmb.svg
	inkscape -w 600 -b white -o $@ $<
	python3 optimize.py $@ $@

ebmb_square.png: ebmb.png
	convert $< -gravity center -extent 600x600 $@

ebmb_banner.png: ebmb.png
	convert $< -gravity center -extent 1280x640 $@

clean:
	rm -f *.aux *.log *.pdf *_*.png

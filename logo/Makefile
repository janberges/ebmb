.PHONY: all clean

w = 15
in = 2.54

all: ebmb.png ebmb_square.png ebmb_banner.png

%.pdf: %.tex
	pdflatex --interaction=batchmode $*

ebmb.png: ebmb.pdf
	convert -density `perl -e "print 640 / ($(w) / $(in))"` $< \
        -flatten PNG8:$@

ebmb_square.png: ebmb.pdf
	convert -density `perl -e "print 640 / ($(w) / $(in))"` $< \
        -gravity center -extent 640x640 PNG8:$@

ebmb_banner.png: ebmb.pdf
	convert -density `perl -e "print 1280 / ($(w) / $(in))"` $< \
        -gravity center -extent 1280x640 PNG8:$@

clean:
	rm -f *.aux *.log *.pdf *.png

#!/bin/bash

w=15
h=6
in=2.54

pdflatex --interaction=batchmode ebmb

convert -density `perl -e "print 640 / ($w / $in)"` ebmb.pdf \
    -flatten PNG8:ebmb.png

convert -density `perl -e "print 640 / ($w / $in)"` ebmb.pdf \
    -background white -gravity center -extent 640x640 PNG8:ebmb_square.png

convert -density `perl -e "print 1280 / ($w / $in)"` ebmb.pdf \
    -background white -gravity center -extent 1280x640 PNG8:ebmb_banner.png

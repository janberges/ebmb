.PHONY: all clean

all: $(patsubst %.py, %.png, $(filter-out a2f.py, $(wildcard *.py)))

clean:
	rm -f *.png

%.png: %.py
	python3 $^

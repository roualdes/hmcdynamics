## Copyright (c) 2019, Edward A. Roualdes
## Distributed under the terms of the Modified BSD License.

all: html

MAINFILE:=report
RMD = $(MAINFILE:=.Rmd)
HTML = $(MAINFILE:=.html)

.PHONY: html clean
html: $(HTML)

%.html: %.Rmd
	Rscript -e "rmarkdown::render('$<')"

clean:
	-rm -f $(HTML)

open:
	open $(HTML)

copy:
	cpy $(CURDIR)/$(HTML)

doc:
	Rscript -e 'devtools::document()'

test:
	Rscript -e 'devtools::test()'

build:
	R CMD build .

check:
	Rscript -e 'devtools::check()'

install:
	Rscript -e 'devtools::install()'

clean:
	rm -rf *.tar.gz *.Rcheck

all:  docs install

docs:
	R -e "devtools::document()"
build:
	(cd ..; R CMD build --no-build-vignettes TrenaProjectAD-PRN)

install:
	(cd ..; R CMD INSTALL --no-test-load TrenaProjectAD-PRN)

check:
	(cd ..; R CMD check `ls -t TrenaProjectAD-PRN) | head -1`)

test:
	for x in inst/unitTests/test_*.R; do echo $$x; R -f $$x; done


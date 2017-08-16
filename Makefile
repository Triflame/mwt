
distribute: clean
	cd ..; tar zcvf mwt.tar.gz mwt
	cp ~/mwt.tar.gz ~/../www/htdocs/tomo07/07_ann/reproducible_codes/chaiwoot/mwt

clean:
	cd data; make clean


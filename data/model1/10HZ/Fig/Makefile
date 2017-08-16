convert: *.vpl
	@for file in *.vpl;do vpconvert format=jpg $$file;done
gif:
	convert -delay 30 -loop 0 *.jpg csg.gif

view:
	for file in *.vpl; do sfpen $$file; done

clean:
	rm -f *.jpg *.gif

.PHONY:
	convert gif view clean

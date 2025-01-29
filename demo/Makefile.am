
# demo Makefile

noinst_SCRIPTS = \
run \
bif.fp \
contour.fp \
controls.fp \
dlm.fp \
gaf.fp \
goland.fp \
latent.fp \
lti1.fp \
nacpv.fp \
nl1.fp \
nl2.fp \
nl3.fp \
nl4.fp \
nl5.fp \
reduction.fp \
reed.fp \
ss.fp \
vacm.fp \
whirl.fp \
whirlnl.fp

EXTRA_DIST = \
$(noinst_SCRIPTS) \
contour.op4 \
contour.uf \
lco.op4 \
lco.uf \
ha145b.op4

clean-local:
		find . -name "valgrind*" -exec rm -f {} \;
		find . -name "*.gprof" -exec rm -f {} \;

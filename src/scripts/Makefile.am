
# src/scripts Makefile

FORCE:

bin_SCRIPTS = \
addheader \
biggest \
bootstrap \
cgrep \
changeversion \
checkVersion \
cleanout \
diffmm \
dog \
fclean \
fig2pdf \
findlink \
fkill \
flapsLog \
flapstest.sh \
fold.sh \
gyroeig.m \
hidden \
hogs \
kind \
ldpackage \
lines \
linked \
lspath \
memmapfile.m \
mgrep \
mklinks \
mminfo.m \
mmread.m \
mmwrite.m \
mmdiff.m \
motd \
nodot \
optfcn.m \
perffix \
rename \
rmtrace \
spy \
targets \
tbig \
tgrep \
traceback \
when \
wraptrace

#EXTRA_DIST = $(bin_SCRIPTS) $(noinst_SCRIPTS)
EXTRA_DIST = $(bin_SCRIPTS)


install-exec-hook:
	touch $(FROOT)/log
	chmod 666 $(FROOT)/log


# doc/figures Makefile
# NOTES:
# - xfig MUST be run in the repo directory, not the build directory,
#   because an xfig save destroys the link and creates the fig file locally
#   
# - figures produced in xfig may use Latex math, e.g.
#   % 12 pt font ---v
#     4 0 0 50 -1 0 12 0.0000 2 150 495 3092 1588 $\\int$\001
#   % Notice the 2 here ------^ which means latex fonts; also the double backslash

# - Figures from xfig are converted to eps with fig2eps

# - To extract pages from manual.pdf and create a new pdf file:
#      texexec --pdfselect --selection=78,79,80 manual.pdf

SUFFIXES = .fig .tex .eps .ps .pdf

FORCE:

# obsolete:
# contour_screen.png
# ffwp.pdf
# goland-contours.png

nosrc = \
controls_vso.png \
goland-contour.jpg \
goland-vso.jpg \
lco_screen.png \
matview.pdf \
mm_matlab.pdf \
viz.png \
viz-amviz.png

# obsolete:
# contour_struct.fig
# controls-vso22-23.fig
# coords.fig
# freeplay.fig
# gaf1.fig
# gaf2.fig
# gapGamma.fig
# ha145b.fig

figs = \
beta.fig \
bif.fig \
bilinear-force.fig \
bma_modes.fig \
bmodes.fig \
cms.fig \
cms_modes.fig \
contour_start.fig \
df.fig \
dlm.fig \
gap-force.fig \
globalcs.fig \
goland.fig \
goland-vso.fig \
gyro.fig \
lti-block.fig \
lti-example.fig \
lti-schematic.fig \
soe.fig \
ss-cms.fig \
ss-km.fig \
ss_phi.fig \
viz-manip.fig \
whirl-vso.fig \
whirl-gaf.fig \
whirl-lbar.fig

# XXX eliminate eps? redundant
eps = \
beta.eps \
bif.eps \
bilinear-force.eps \
bma_modes.eps \
bmodes.eps \
cms.eps \
cms_modes.eps \
contour_start.eps \
df.eps \
dlm.eps \
gap-force.eps \
globalcs.eps \
goland.eps \
goland-vso.eps \
gyro.eps \
lti-block.eps \
lti-example.eps \
lti-schematic.eps \
soe.eps \
ss-cms.eps \
ss-km.eps \
ss_phi.eps \
viz-manip.eps \
vso2d/vso2d.eps \
whirl-vso.eps \
whirl-gaf.eps \
whirl-lbar.eps

# XXX need to list *only* the .fig files, then make a list of .eps files from that

%.eps:	%.fig
	fig2eps --packages=amsmath,graphicx,color --nogv $<

all:	$(eps)

# list of stuff to include in dist tarfile XXX no eps?
EXTRA_DIST = $(nosrc) $(figs) $(eps)

# noinst_DATA = \
# $(figs) \
# $(texfigs)

list:		FORCE
	sed -n 's/^.*includegraphics.*{\([^}]*\)}.*$/\1/' ../manual.tex

tarfile:		FORCE
	tar zcf manual.tgz manual.tex fig *.eps Makefile


CLEANFILES = \
$(eps)

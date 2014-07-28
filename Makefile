
OX := "C:\Program Files (x86)\OxMetrics6\Ox\bin\oxl.exe"
OXFLAGS := 
COPY := copy
OXDOC := "C:\Program Files (x86)\oxdoc\bin\oxdoc.bat"
ERASE := erase

vpath %.ox .
vpath %.h .
vpath %.ox.html .

.PHONY : document
document:
	$(OXDOC) *.ox
	
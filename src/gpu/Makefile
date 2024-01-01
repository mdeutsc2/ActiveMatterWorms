# Makefile for compiling amatter2.chpl with Chapel compiler

CHPL = chpl
SRCDIR = src/
#SOURCE = $(SRCDIR)amatter.chpl
#SOURCE = $(SRCDIR)amatter3d.chpl
SOURCE = amatter_gpu.chpl
BUILDDIR = build/
OUTPUT = amatter.x

GPUFLAGS = --fast --optimize --vectorize --inline --print-callstack-on-error --gpu-specialization --ccflags=--cuda-path=/usr/lib/cuda --report-gpu

#GPUFLAGS = --debug --print-callstack-on-error --gpu-specialization --ccflags=--cuda-path=/usr/lib/cuda --report-gpu
clean:
	rm -f $(OUTPUT) *.x

cuda:
	$(CHPL) $(GPUFLAGS) $(SOURCE) -o $(OUTPUT)

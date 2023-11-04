# Makefile for compiling amatter2.chpl with Chapel compiler

CHPL = chpl
SRCDIR = src/
SOURCE = $(SRCDIR)amatter.chpl
#SOURCE = $(SRCDIR)amatter_new_gpu.chpl
BUILDDIR = build/
OUTPUT = amatter.x

CHPLFLAGS = --fast --savec=$(BUILDDIR) --optimize --vectorize --specialize --optimize-loop-iterators --print-passes --print-commands --print-callstack-on-error
DEBUGFLAGS = -g --savec=$(BUILDDIR) --ldflags -v --print-passes --print-commands --checks --print-callstack-on-error
GPUFLAGS = --fast --optimize --ldflags -v --print-passes --print-commands --gpu-specialization --ccflags=--cuda-path=/usr/lib/cuda --report-gpu

all: $(OUTPUT)

$(OUTPUT): $(SOURCE)
	$(CHPL) $(CHPLFLAGS) $(SOURCE) -o $(OUTPUT)

clean:
	rm -f $(OUTPUT) *.x

debug: $(SOURCE)
	$(CHPL) $(DEBUGFLAGS) $(SOURCE) -o $(OUTPUT)

cuda:
	$(CHPL) $(GPUFLAGS) $(SOURCE) -o $(OUTPUT)

run:
	set -e; \
	dirname=data_$$(date +%d_%m_%Y); \
	if [ -d "$$dirname" ]; then \
		count=2; \
		while [ -d "$$dirname-$$count" ]; do \
			count=$$(($$count+1)); \
		done; \
		dirname="$$dirname-$$count"; \
	fi; \
	mkdir -p $$dirname; \
	cd $$dirname && time ../$(OUTPUT); \
	output_file=amatter_$$(date +%d_%m_%Y); \
	if [ -n "$$count" ]; then \
		output_file=$$output_file-$$count; \
	fi; \
	output_file=$$output_file.xyz; \
	python3 ../merge_xyz.py $$output_file; \
	cp $$output_file ..

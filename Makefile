# Makefile for compiling amatter2.chpl with Chapel compiler

CHPL = chpl
CHPLFLAGS = --fast --optimize --ldflags -v --print-passes --print-commands --print-callstack-on-error
DEBUGFLAGS = -g --savec=c_srcs --ldflags -v --print-passes --print-commands --checks --print-callstack-on-error
GPUFLAGS = --fast --optimize --ldflags -v--print-passes --print-commands --gpu-arch=sm60 --gpu-specialization

# SRCS = amatter.chpl
# VPATH = src/
SRCDIR = src/
SOURCE = $(SRCDIR)amatter.chpl
# SOURCE := $(SRCS)
OUTPUT = amatter.x
SETUP_SCRIPT = /mnt/SCRATCH/chapel-1.31.0/util/setchplenv.bash
all: $(OUTPUT)

$(OUTPUT): $(SOURCE)
	$(CHPL) $(CHPLFLAGS) $(SOURCE) -o $(OUTPUT)

clean:
	rm -f $(OUTPUT) *.x

debug: $(SOURCE)
	$(CHPL) $(DEBUGFLAGS) $(SOURCE) -o $(OUTPUT)

setup:
	source $(SETUP_SCRIPT)

run:
	set -e; \
	if [ ! -e "$(OUTPUT)" ]; then \
		$(MAKE) all; \
	fi; \
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
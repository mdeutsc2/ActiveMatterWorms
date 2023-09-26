# Makefile for compiling simple.chpl with Chapel compiler

CHPL = chpl
CHPLFLAGS = --fast --optimize --print-passes --print-commands
DEBUGFLAGS = -g --savec=c_srcs --print-passes --print-commands --checks --print-callstack-on-error
GPUFLAGS = --fast --optimize --print-passes --print-commands --gpu-arch=sm60 --gpu-specialization

SOURCE = simple.chpl
OUTPUT = simple.x
SETUP_SCRIPT = /mnt/SCRATCH/chapel-1.31.0/util/setchplenv.bash
all: $(OUTPUT)

$(OUTPUT): $(SOURCE)
	$(CHPL) $(CHPLFLAGS) $(SOURCE) -o $(OUTPUT)

debug: $(SOURCE)
	$(CHPL) $(DEBUGFLAGS) $(SOURCE) -o $(OUTPUT)

cuda: $(SOURCE)
	CHPL_GPU=nvidia $(CHPL) $(GPUFLAGS) $(SOURCE) -o $(OUTPUT)
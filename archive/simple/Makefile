# Makefile for compiling simple.chpl with Chapel compiler

CHPL = chpl
CHPLFLAGS = --fast --optimize --ldflags -v --print-passes --print-commands --print-callstack-on-error
DEBUGFLAGS = -g --savec=c_srcs --ldflags -v --print-passes --print-commands --checks --print-callstack-on-error
GPUFLAGS = --fast --optimize --ldflags -v--print-passes --print-commands --gpu-arch=sm60 --gpu-specialization

SOURCE = simple.chpl
OUTPUT = simple.x
SETUP_SCRIPT = /mnt/SCRATCH/chapel-1.31.0/util/setchplenv.bash
all: $(OUTPUT)

$(OUTPUT): $(SOURCE)
	$(CHPL) $(CHPLFLAGS) $(SOURCE) -o $(OUTPUT)

debug: $(SOURCE)
	$(CHPL) $(DEBUGFLAGS) $(SOURCE) -o $(OUTPUT)

clean:
	rm *.x

cuda: simple2.chpl
	CHPL_GPU=nvidia $(CHPL) $(CHPLFLAGS) simple2.chpl -o simple2.x

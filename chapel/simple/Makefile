# Makefile for compiling simple.chpl with Chapel compiler

CHPL = chpl
CHPLFLAGS = --fast --optimize --print-passes --print-commands

SOURCE = simple.chpl
OUTPUT = simple.x
SETUP_SCRIPT = /mnt/SCRATCH/chapel-1.31.0/util/setchplenv.bash
all: $(OUTPUT)

$(OUTPUT): $(SOURCE)
	$(CHPL) $(CHPLFLAGS) $(SOURCE) -o $(OUTPUT)

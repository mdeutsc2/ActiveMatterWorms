OUTPUT = main.x
SOURCE = main.chpl
CHPL = chpl

CHPLFLAGS = --print-callstack-on-error --ccflags=-fopenmp=libomp

all: $(OUTPUT)

$(OUTPUT): $(SOURCE)
	$(CHPL) $(CHPLFLAGS) $(SOURCE) -o $(OUTPUT)


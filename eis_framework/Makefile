# EIS Framework Makefile
# Independent from main LMFIT project

CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c99 -Isrc -I../lib
LDFLAGS = -L../lib -llmmin -lm

# Directories
SRCDIR = src
MODELDIR = models
EXAMPLEDIR = examples
TESTDIR = tests
BUILDDIR = build

# Source files
CORE_SOURCES = $(SRCDIR)/eis_core.c $(SRCDIR)/eis_fitting.c $(SRCDIR)/eis_output.c
MODEL_SOURCES = $(wildcard $(MODELDIR)/*.c)
EXAMPLE_SOURCES = $(wildcard $(EXAMPLEDIR)/*.c)
TEST_SOURCES = $(wildcard $(TESTDIR)/*.c)

# Object files
CORE_OBJECTS = $(CORE_SOURCES:%.c=$(BUILDDIR)/%.o)
MODEL_OBJECTS = $(MODEL_SOURCES:%.c=$(BUILDDIR)/%.o)
EXAMPLE_OBJECTS = $(EXAMPLE_SOURCES:%.c=$(BUILDDIR)/%.o)
TEST_OBJECTS = $(TEST_SOURCES:%.c=$(BUILDDIR)/%.o)

# Library and executables
LIBRARY = $(BUILDDIR)/libeis.a
DEMO_EXEC = $(BUILDDIR)/eis_demo
TEST_EXEC = $(BUILDDIR)/eis_test

# Create build directory
$(shell mkdir -p $(BUILDDIR)/$(SRCDIR) $(BUILDDIR)/$(MODELDIR) $(BUILDDIR)/$(EXAMPLEDIR) $(BUILDDIR)/$(TESTDIR))

.PHONY: all library demo test clean install help

all: library demo

# Build static library
library: $(LIBRARY)

$(LIBRARY): $(CORE_OBJECTS) $(MODEL_OBJECTS)
	ar rcs $@ $^
	@echo "EIS library built: $(LIBRARY)"

# Build demonstration program
demo: $(DEMO_EXEC)

$(DEMO_EXEC): $(EXAMPLE_OBJECTS) $(LIBRARY)
	$(CC) $(EXAMPLE_OBJECTS) -L$(BUILDDIR) -leis $(LDFLAGS) -o $@
	@echo "Demo executable built: $(DEMO_EXEC)"

# Build test suite
test: $(TEST_EXEC)

$(TEST_EXEC): $(TEST_OBJECTS) $(LIBRARY)
	$(CC) $(TEST_OBJECTS) -L$(BUILDDIR) -leis $(LDFLAGS) -o $@
	@echo "Test executable built: $(TEST_EXEC)"

# Pattern rules for object files
$(BUILDDIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Dependencies
$(CORE_OBJECTS) $(MODEL_OBJECTS) $(EXAMPLE_OBJECTS): $(SRCDIR)/eis_framework.h

# Run tests
run-test: test
	@echo "Running EIS framework tests..."
	./$(TEST_EXEC)

# Run demo
run-demo: demo
	@echo "Running EIS framework demonstration..."
	./$(DEMO_EXEC)

# Run with sample data
run-sample: demo
	@echo "Running demo with sample data..."
	./$(DEMO_EXEC) examples/sample_data.txt results.txt

# Clean build files
clean:
	rm -rf $(BUILDDIR)
	rm -f *.txt *.dat

# Install library and headers
install: library
	@echo "Installing EIS framework..."
	sudo cp $(LIBRARY) /usr/local/lib/
	sudo cp $(SRCDIR)/eis_framework.h /usr/local/include/
	sudo ldconfig
	@echo "Installation complete"

# Uninstall
uninstall:
	sudo rm -f /usr/local/lib/libeis.a
	sudo rm -f /usr/local/include/eis_framework.h
	sudo ldconfig

# Help
help:
	@echo "EIS Framework Build System"
	@echo "=========================="
	@echo ""
	@echo "Targets:"
	@echo "  all        - Build library and demo (default)"
	@echo "  library    - Build EIS static library"
	@echo "  demo       - Build demonstration program"
	@echo "  test       - Build test suite"
	@echo "  run-demo   - Build and run demonstration"
	@echo "  run-test   - Build and run tests"
	@echo "  run-sample - Run demo with sample data"
	@echo "  clean      - Remove all build files"
	@echo "  install    - Install library system-wide"
	@echo "  uninstall  - Remove installed files"
	@echo "  help       - Show this help"
	@echo ""
	@echo "Usage examples:"
	@echo "  make                     # Build everything"
	@echo "  make run-demo            # Quick test run"
	@echo "  make run-sample          # Test with sample data"
	@echo "  ./build/eis_demo my_data.txt results.txt"
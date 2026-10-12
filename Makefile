# --- Compiler selection (default: g++) ---
# Use `make CXX=clang++` to switch to clang++.
CXX ?= g++

CXXFLAGS ?= -std=c++20 -Wall -Wextra -O3 -DNDEBUG -mtune=native -march=native -funroll-loops
CPPFLAGS += -Iinclude

ifeq ($(BIG_N),1)
	CPPFLAGS += -DUSE_128BIT_PACKING
MODE = large
else
MODE = standard
endif

ifeq ($(SCALAR),1)
CPPFLAGS += -DKLC_DISABLE_SIMD
MODE := $(MODE)-scalar
endif

# --- Variables ---
SRCS = main.cpp src/Polynomial.cpp src/Permutation.cpp src/Bruhat.cpp src/KLCalculator.cpp src/Lascoux.cpp src/Commands.cpp src/Cli.cpp
OBJDIR = obj/$(MODE)
OBJS = $(addprefix $(OBJDIR)/,$(notdir $(SRCS:.cpp=.o)))
TARGET ?= klc_cpp
TEST_TARGET = $(OBJDIR)/core_tests
CORE_OBJS = $(filter-out $(OBJDIR)/main.o $(OBJDIR)/Cli.o $(OBJDIR)/Commands.o,$(OBJS))

.PHONY: all clean test

# --- Main Build ---
all: $(TARGET)

$(TARGET): $(OBJS) FORCE
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

INCLUDES = $(wildcard include/*.hpp include/detail/*.hpp)

$(OBJDIR)/%.o: src/%.cpp $(INCLUDES) | $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(OBJDIR)/main.o: main.cpp $(INCLUDES) | $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(TEST_TARGET): tests/core_tests.cpp $(CORE_OBJS) $(INCLUDES)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(CORE_OBJS)

test: $(TARGET) $(TEST_TARGET)
	./$(TEST_TARGET)
	./$(TARGET) --verify verify.txt
	python tests/cli_tests.py ./$(TARGET)
	python tests/tools_tests.py

.PHONY: FORCE
FORCE:

-include $(OBJS:.o=.d)

# --- Cleanup ---
clean:
	rm -rf obj $(TARGET) $(TARGET).exe

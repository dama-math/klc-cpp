# --- Compiler selection (default: clang++) ---
# Use `make COMPILER=gcc` or `make CXX=g++` to switch to g++.
COMPILER ?=
ifeq ($(COMPILER),gcc)
	CXX := g++
else
	CXX := clang++
endif

CXXFLAGS = -std=c++20 -Wall -O3 -mtune=native -march=native -funroll-loops -Iinclude

ifeq ($(BIG_N),1)
	CXXFLAGS += -DUSE_128BIT_PACKING
endif

# --- Variables ---
SRCS = main.cpp src/Polynomial.cpp src/Permutation.cpp src/Bruhat.cpp src/KLCalculator.cpp src/Lascoux.cpp src/Commands.cpp src/Utils.cpp
OBJDIR = obj
OBJS = $(addprefix $(OBJDIR)/,$(notdir $(SRCS:.cpp=.o)))
TARGET = klc_cpp

.PHONY: all clean

# --- Main Build ---
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OBJDIR)/%.o: src/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR)/main.o: main.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

# --- Cleanup ---
clean:
	rm -rf $(OBJDIR) $(TARGET)

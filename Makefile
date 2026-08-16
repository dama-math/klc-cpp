# --- Compiler selection (default: g++) ---
# Use `make CXX=clang++` to switch to clang++.
CXX ?= g++

CXXFLAGS = -std=c++20 -Wall -O3 -mtune=native -march=native -funroll-loops -Iinclude

ifeq ($(BIG_N),1)
	CXXFLAGS += -DUSE_128BIT_PACKING
endif

# --- Variables ---
SRCS = main.cpp src/Polynomial.cpp src/Permutation.cpp src/Bruhat.cpp src/KLCalculator.cpp src/Lascoux.cpp src/Commands.cpp
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

CXX=g++
CXXFLAGS=-g -std=c++17 -O1 -fopenmp  
#-fsanitize=address 

EXE=render
ODIR=build

SRCS = $(wildcard *.cpp) $(wildcard */*.c)
DEPS = $(wildcard *.h) 

render: $(SRCS) $(DEPS)
	$(CXX) -o render $(SRCS) $(CXXFLAGS) $(LIBS)

.PHONY: incremental
incremental: $(OBJS)
	$(CXX) -o render $^ $(CXXFLAGS) $(LIBS)

.PHONY: clean
clean:
	rm -r $(OBJS) $(EXE)

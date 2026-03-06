#CXX = g++
CXX = clang++
CXXFLAGS = -std=c++17 -fopenmp -I/usr/include/eigen3 -Wall -g -MMD -MP
TARGET = ns_fda_3D
SRCS = src/main.cpp src/inout.cpp src/mesh_n_model.cpp src/compute_flow.cpp src/generated.cpp
OBJS = $(SRCS:.cpp=.o)

.PHONY: clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)
    
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

-include $(OBJS:.o=.d)

clean:
	rm -f $(TARGET) $(OBJS) $(OBJS:.o=.d)

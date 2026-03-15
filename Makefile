# Makefile 
CXX = g++
CXXFLAGS = -std=c++17 -O7

SRCS = run.cpp RiemannSolvers.cpp IO.cpp DragIntegrators.cpp BoundaryConditions.cpp Reconstruction.cpp
OBJS = run.o RiemannSolvers.o IO.o DragIntegrators.o BoundaryConditions.o Reconstruction.o 
HEADERS = RiemannSolvers.h IO.h classes.h DragIntegrators.h BoundaryConditions.h Reconstruction.h indices.h

TARGET = pigpen.exe

# Default target
all: $(TARGET)

# Link object files to create the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Pattern rule for compiling .cpp files to .o files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean target
clean:
	rm -f $(TARGET) $(OBJS)
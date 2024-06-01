# Compiler and flags
NVCC = nvcc
GCC = gcc
CFLAGS =

# Source files
CU_SRCS = dev_main_cuda.cu
C_SRCS = cuda_utils.c

# Object files
CU_OBJS = $(CU_SRCS:.cu=.o)
C_OBJS = $(C_SRCS:.c=.o)

# Executable
TARGET = dev_main_cuda

# Default rule
all: $(TARGET)

# Rule to link the final executable
$(TARGET): $(CU_OBJS) $(C_OBJS)
	$(NVCC) -o $@ $(CU_OBJS) $(C_OBJS)

# Rule to compile CUDA source files
%.o: %.cu
	$(NVCC) -c -o $@ $<

# Rule to compile C source files
%.o: %.c
	$(GCC) $(CFLAGS) -c -o $@ $<

# Clean rule
clean:
	rm -f $(CU_OBJS) $(C_OBJS) $(TARGET)

# Phony targets
.PHONY: all clean

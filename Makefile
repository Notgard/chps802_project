#
# MAIN CONFIGURATION (to configure)
#

EXEC = main gen_matrix stats
PROFILE_EXEC=main
OBJECTS = utils.o omp_utils.o matmul.o
PROJECT_NAME = gaussian_project
PROFILE_FILE=generated_matrix.txt

#
# SUFFIXES (must not change it)
#

.SUFFIXES: .c .o

#
# OBJECTS (must not change it)
#

EXEC_O = $(EXEC:=.o)
OBJECTS_O = $(OBJECTS) $(EXEC_O)

#
# ARGUMENTS AND COMPILER (to configure)
#

CC = gcc
OPTIMIZER_FLAGS = -O2 -pg -g 
CCFLAGS_STD = -fopenmp -pg -g -Wall -Wextra -Wshadow #$(OPTIMIZER_FLAGS)
CCFLAGS_DEBUG = -D _DEBUG_
CCFLAGS_OMP = -D _OMP_
CCFLAGS = $(CCFLAGS_STD) -pg -g 
CCLIBS = -lm -fopenmp -pg -g 

#
# RULES (must not change it)
#

all: msg $(OBJECTS) $(EXEC_O)
	@echo "Create executables..."
	@for i in $(EXEC); do \
	$(CC) -o $$i $$i.o $(OBJECTS) $(CCLIBS); \
	done
	@echo "Done."

msg:
	@echo "Create objects..."

debug: CCFLAGS = $(CCFLAGS_STD) $(CCFLAGS_DEBUG)
debug: all

omp: CCFLAGS = $(CCFLAGS_STD) $(CCFLAGS_OMP)
omp: all

#
# DEFAULT RULES (must not change it)
#

%.o : %.c
	@cd $(dir $<) && ${CC} ${CCFLAGS} -c $(notdir $<) -o $(notdir $@)

#
# MAIN RULES (must not change it)
#

# You can add your own commands
clean:
	@echo "Delete objects, temporary files..."
	@rm -f $(OBJECTS) $(EXEC_O)
	@rm -f *~ *#
	@rm -f $(EXEC)
	@rm -f dependancies
	@rm -f *gmon.out
	@echo "Done."

# DEPENDANCIES
# This section is completed automatically

# Profiling rule
profile: all
	@echo "Running executables for profiling..."
	@for i in $(PROFILE_EXEC); do \
		echo "Profiling $$i..."; \
		./$$i $(PROFILE_FILE); \
		if [ -f gmon.out ]; then \
			gprof $$i gmon.out > $$i.prof; \
			mv gmon.out $$i.gmon.out; \
		else \
			echo "gmon.out not found for $$i."; \
		fi \
	done
	@echo "Profiling done."


##################################################################
# Import machine specific compilers, options, flags, etc.. 
##################################################################

BRAID_DIR=../braid
include ../makefile.inc

##################################################################
# Build examples 
##################################################################

BRAID_FLAGS = -I$(BRAID_DIR)
BRAID_LIB_FILE = $(BRAID_DIR)/libbraid.a

C_NOHYPRE = xbraid-SWE

.PHONY: all clean cleanout

.SUFFIXES:
.SUFFIXES: .c .cpp

# put this rule first so it becomes the default
all: $(C_NOHYPRE) $(C_EXAMPLES) $(CXX_EXAMPLES)

# Rule for compiling .c files
%: %.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) -fPIC $(CFLAGS) $(BRAID_FLAGS) $(@).c -o $@\
 	$(BRAID_LIB_FILE) $(LFLAGS)

# Rule for compiling .cpp files; links with mfem
%: %.cpp $(BRAID_LIB_FILE) 
	@echo "Building" $@ "..."
	$(MPICXX) -fPIC $(CXXFLAGS) $(BRAID_FLAGS) \
	$< -o $@ $(BRAID_LIB_FILE) $(LFLAGS)

clean: cleanout
	rm -f *.o $(C_NOHYPRE) $(F_NOHYPRE) $(C_EXAMPLES) $(CXX_EXAMPLES) $(F_EXAMPLES) *ror_norm* *_err_* *_mesh* *_sol_*
	rm -rf *.dSYM

cleanout:
	rm -f ex*.out.*


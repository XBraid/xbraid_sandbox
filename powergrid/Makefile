MPICC=mpicc

BRAID_DIR=$(HOME)/Numerics/xbraid/braid
BRAID_FLAGS = -I$(BRAID_DIR)
BRAID_LIB_FILE = $(BRAID_DIR)/libbraid.a

LFLAGS=-lm

C_NOHYPRE = ex-01-expanded-powergrid ex-01-expanded-powergrid-aAsSeperateVar ex-01-expanded-powergrid-aAsSeperateVar2 ex-01-expanded-powergrid-doublevector ex-01-expanded-powergrid-optim

.PHONY: all clean

.SUFFIXES:
.SUFFIXES: .c

# put this rule first so it becomes the default
all: $(C_NOHYPRE) 

# Rule for building ex-01-expanded-powergrid
ex-01-expanded-powergrid: ex-01-expanded-powergrid.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) $(@).c -o $@ $(BRAID_LIB_FILE) $(LFLAGS)

# Rule for building ex-01-expanded-powergrid
ex-01-expanded-powergrid-doublevector: ex-01-expanded-powergrid-doublevector.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) $(@).c -o $@ $(BRAID_LIB_FILE) $(LFLAGS)

# Rule for building ex-01-expanded-powergrid
ex-01-expanded-powergrid-aAsSeperateVar: ex-01-expanded-powergrid-aAsSeperateVar.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) $(@).c -o $@ $(BRAID_LIB_FILE) $(LFLAGS)


# Rule for building ex-01-expanded-powergrid
ex-01-expanded-powergrid-aAsSeperateVar2: ex-01-expanded-powergrid-aAsSeperateVar2.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) $(@).c -o $@ $(BRAID_LIB_FILE) $(LFLAGS)
#
# Rule for building ex-01-expanded-powergrid-optim
ex-01-expanded-powergrid-optim: ex-01-expanded-powergrid-optim.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) $(@).c -o $@ $(BRAID_LIB_FILE) $(LFLAGS)


# Rule for compiling .c files
%: %.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) $(HYPRE_FLAGS) $(@).c -o $@\
 $(BRAID_LIB_FILE) $(HYPRE_LIB) $(LFLAGS)


clean: 
	rm -f *.o $(C_NOHYPRE) *ror_norm* *_err_* *_mesh* *_sol_*
	rm -rf *.dSYM


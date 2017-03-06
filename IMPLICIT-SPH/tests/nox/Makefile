info :
	@echo "make one TEST=scalar"

include ../Makefile.bread

TEST = scalar_aztecOO
TESTS = scalar_aztecOO

OBJ = 

LDFLAGS = 

CPP = g++
CFLAGS = -g -O0


CC_WORK = $(CC) $(CCFLAGS) $(INCLUDES)


one : $(OBJ) test_$(TEST).o
	@echo "Linking $<" 
	$(CC_WORK) -o $(TEST) test_$(TEST).o $(OBJ) $(LIBS) $(LDFLAGS)

%.o : %.cpp
	@echo "Compiling $<" 
	$(CC_WORK) -o $@ -c $<

clean: 
	@echo $(CC_WORK)
	/bin/rm -f *.o *~ 
	for t in $(TESTS) ; do rm -f $$t ; done

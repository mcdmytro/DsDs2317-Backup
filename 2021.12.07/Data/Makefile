override LINKFLAGS += -shared -L$(BELLE_RUN_DIR)/lib/so -L$(BELLE_RUN_DIR)/bin

override DEBUG =-g -Wall
override COMPILEFLAGS += -I$(BELLE_TOP_DIR)/include -DHEP_SHORT_NAMES  -fpic -pipe
SRC = reco2460.cc
# If there are several .cc files to be compiled you must point them above

%.o: %.cc %.h
	gcc -c $(COMPILEFLAGS) $(DEBUG) -MMD $(INC) -o $@ $<

User_reco.so: reco2460.o myutl.o
	gcc $(LINKFLAGS) $(DEBUG)    $^ -o $@
clean:
	-rm *.o *.so *.d


-include $(DEP)

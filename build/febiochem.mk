
SRC = $(wildcard $(FCDIR)FEBioChem/*.cpp)
OBJ = $(patsubst $(FCDIR)FEBioChem/%.cpp, %.o, $(SRC))
DEP = $(patsubst $(FCDIR)FEBioChem/%.cpp, %.d, $(SRC))


SO = libfebiochem_$(PLAT).$(SFX)
LIB = $(FCDIR)build/lib/$(SO)

FECORE = $(FEBLIB)/libfecore_$(PLAT).a

FEBIOMECH = $(FEBLIB)/libfebiomech_$(PLAT).a

FEBIOLIBS = $(FEBIOMECH) $(FECORE)

$(LIB): $(OBJ)
ifeq ($(findstring lnx,$(PLAT)),lnx)
		$(CC) $(LNKFLG) -shared -Wl,-soname,$(SO) -o $(LIB) $(OBJ) $(FEBIOLIBS) $(INTEL_LIB)libiomp5.a
else ifeq ($(findstring gcc,$(PLAT)),gcc)
		$(CC) $(LNKFLG) -shared -Wl,-soname,$(SO) -o $(LIB) $(OBJ) $(FEBIOLIBS) $(INTEL_LIB)libiomp5.a
else
		$(CC) -dynamiclib $(FLG) -o $(LIB) $(OBJ) $(FEBIOLIBS) $(INTEL_LIB)libiomp5.a
endif

%.o: $(FCDIR)FEBioChem/%.cpp
	$(CC) $(INC) $(DEF) $(FLG) -MMD -c -o $@ $<

clean:
	$(RM) *.o *.d $(LIB)

-include $(DEP)

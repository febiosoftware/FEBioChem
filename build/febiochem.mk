
SRC = $(wildcard $(FCDIR)source/*.cpp)
OBJ = $(patsubst $(FCDIR)source/%.cpp, %.o, $(SRC))
DEP = $(patsubst $(FCDIR)source/%.cpp, %.d, $(SRC))


SO = libfebiochem_$(PLAT).$(SFX)
LIB = $(FCDIR)build/lib/$(SO)

FECORE = $(FEBLIB)/libfecore_$(PLAT).a

FEBIOMECH = $(FEBLIB)/libfebiomech_$(PLAT).a

FEBIOLIBS = $(FEBIOMECH) $(FECORE)

$(LIB): $(OBJ)
ifeq ($(findstring lnx,$(PLAT)),lnx)
		$(CC) $(LNKFLG) -shared -Wl,-soname,$(SO) -o $(LIB) $(OBJ) $(FEBIOLIBS)
else ifeq ($(findstring gcc,$(PLAT)),gcc)
		$(CC) $(LNKFLG) -shared -Wl,-soname,$(SO) -o $(LIB) $(OBJ) $(FEBIOLIBS)
else
		$(CC) -dynamiclib $(FLG) -o $(LIB) $(OBJ) $(FEBIOLIBS)
endif

%.o: $(FCDIR)source/%.cpp
	$(CC) $(INC) $(DEF) $(FLG) -MMD -c -o $@ $<

clean:
	$(RM) *.o *.d $(LIB)

-include $(DEP)

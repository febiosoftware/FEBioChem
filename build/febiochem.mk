THIS = FEBioChem
SRC = $(wildcard $(FCDIR)$(THIS)/*.cpp)
OBJ = $(patsubst $(FCDIR)$(THIS)/%.cpp, %.o, $(SRC))
DEP = $(patsubst $(FCDIR)$(THIS)/%.cpp, %.d, $(SRC))


SO = lib$(PLGN)_$(PLAT).$(SFX)
LIB = $(FCDIR)build/lib/$(SO)

FECORE = $(FEBLIB)/libfecore.a

FEBIOMECH = $(FEBLIB)/libfebiomech.a

FEBIOMIX = $(FEBLIB)/libfebiomix.a

FEBIOLIBS = $(FEBIOMIX) $(FEBIOMECH) $(FECORE)

$(LIB): $(OBJ)
ifeq ($(findstring lnx,$(PLAT)),lnx)
		$(CC) $(LNKFLG) -shared -Wl,-soname,$(SO) -o $(LIB) $(OBJ) $(FEBIOLIBS) $(INTEL_LIB) -liomp5
else ifeq ($(findstring gcc,$(PLAT)),gcc)
		$(CC) $(LNKFLG) -shared -Wl,-soname,$(SO) -o $(LIB) $(OBJ) $(FEBIOLIBS)
else ifeq ($(findstring osx,$(PLAT)),osx)
		$(CC) -install_name @rpath/lib$(PLGN).dylib -dynamiclib $(FLG) -o $(LIB) $(OBJ) $(FEBIOLIBS) -L/usr/local/opt/libomp/lib -lomp
else
		$(CC) -dynamiclib $(FLG) -o $(LIB) $(OBJ) $(FEBIOLIBS) $(INTEL_LIB) -liomp5
endif

%.o: $(FCDIR)$(THIS)/%.cpp
	$(CC) $(INC) $(DEF) $(FLG) -MMD -c -o $@ $<

clean:
	$(RM) *.o *.d $(LIB)

-include $(DEP)

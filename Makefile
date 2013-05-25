PREFIX  = $(HOME)
TOTAL   = $(HOME)
LIBNAME = NEUS
LIBRARY = lib$(LIBNAME).so

OPT2       = -O2

ROOTCONFIG = root-config

ARCH      := $(shell $(ROOTCONFIG) --arch)
ALTCXX    := $(shell $(ROOTCONFIG) --cxx)
ROOTCFLAGS:= $(shell $(ROOTCONFIG) --cflags)
ROOTLIBS  := $(shell $(ROOTCONFIG) --libs)

ARCHOK     = no

ifeq ($(ARCH),linux)
  CXX      = g++
  CXXFLAGS = $(OPT2) -Wall -fPIC
  SOFLAGS  = -shared
  ARCHOK   = yes
endif

ifeq ($(ARCH),linuxx8664gcc)
  CXX      = g++
  CXXFLAGS = $(OPT2) -Wall -fPIC
  SOFLAGS  = -shared
  ARCHOK   = yes
endif

ifeq ($(ARCH),win32gcc)
  CXX      = g++
  CXXFLAGS = $(OPT2) -Wall
  SOFLAGS  = -shared
  ARCHOK   = yes
endif

# If arch != the above options, a error message is given
ifeq ($(ARCHOK),no)
  $(error $(ARCH) invalid architecture)
endif

# In case that ALTCXX = 0, g++ will be used as CXX
ifneq ($(ALTCXX),)
  CXX = $(ALTCXX)
endif

# Finally, define CXXFLAGS & LIBS
CXXFLAGS+= $(ROOTCFLAGS)
CXXFLAGS+= -g
LIBS     = $(ROOTLIBS)


# Define things related to rootcint
# =================================

ROOTCINT = rootcint

ROOTIFIED_SOURCE := $(LIBNAME)Dict.cc
ROOTIFIED_HEADER := $(ROOTIFIED_SOURCE:.cc=.h)
ROOTIFIED_OBJECT := $(ROOTIFIED_SOURCE:.cc=.o)


# Define SOURCES, HEADERS & OBJECTS 
# ==========================================

SOURCES = $(filter-out $(ROOTIFIED_SOURCE), $(wildcard *.cc))
HEADERS = $(SOURCES:.cc=.h)
OBJECTS = $(SOURCES:.cc=.o)
DEPFILE = $(SOURCES:.cc=.d)

LINKDEF = LinkDef.h


# Define ROOTMAP & variables to create them
# ==================================================

ROOTMAP = $(LIBRARY:.so=.rootmap)

SYMBOLS = `nm -CPu $(LIBRARY) |\
	  awk -F: '/^T/{printf("(%s)\n",$$1)}' |\
	  sort -u | tr '\n' '|'`
ALLLIBS = $(wildcard $(shell $(ROOTCONFIG) --libdir)/*.so)
DEPENDS = $(shell symbols=$(SYMBOLS); \
	  for so in $(ALLLIBS); \
	  do nm -CD --defined-only $$so |\
	  grep -E 'T ('$$symbols')::' > /dev/null &&\
	  echo $$so;\
	  done | sort -u | tr '\n' ' ')

RLIBMAP = rlibmap


# Action starts
# =============

# the first target is the default target, it depends on $(ROOTMAP)
# before "make all", make will include all other makefiles specified
# by the include command
all: $(ROOTMAP)
	@echo
	@echo "* Done!"
	@echo 

# include *.d files, which are makefiles defining dependencies between files
ifeq ($(filter info clean tags,$(MAKECMDGOALS)),)
  -include $(DEPFILE)
endif

# rules to create *.d files
%.d:%.cc
	@echo creating $@
	@set -e; rm -f $@; \
	  $(CXX) -MM $(CPPFLAGS) $< > $@.$$$$; \
	  sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	  rm -f $@.$$$$ 

# lib$(LIBNAME).rootmap can only be created after the creation of lib$(LIBNAME).so. It
# tells ROOT the dependence among libraries. Putting it along with the
# corresponding library allows one to use in CINT the functions defined in the
# library without calling gSystem->Load("lib.so")
$(ROOTMAP): $(LIBRARY)
	@echo
	@echo "* Creating rootmap file:"
	$(RLIBMAP) -o $(ROOTMAP) -l $(LIBRARY) -d $(DEPENDS) -c $(LINKDEF)

# lib$(LIBNAME).so depends on all *.o files.
#  The flag "-shared" is used to create shared libs
#  $@ represents the target, that is, lib$(LIBNAME).so
#  $^ represents all the prerequisites, i.e., all *.o files
$(LIBRARY): $(ROOTIFIED_OBJECT) $(OBJECTS)
	@echo
	@echo "* Creating shared library:"
	$(CXX) $(CXXFLAGS) $(LIBS) $(SOFLAGS) -o $@ $^

# An xxx.o file depends on xxx.cc. It is created with the command:
# 	g++ -c xxx.cc -o xxx.o
# Since this is obvious, "make" automatically does it. 
# There is no need to explicitly write down the rules to do it.

# To use classes & functions directly in ROOT. One has to do the following
$(ROOTIFIED_SOURCE): $(HEADERS) $(LINKDEF)
	@echo 
	@echo "* Rootifying files:" 
	@rm -f $(ROOTIFIED_SOURCE) $(ROOTIFIED_HEADER) 
	$(ROOTCINT) $(ROOTIFIED_SOURCE) -c -p $(CXXFLAGS) $(HEADERS) $(LINKDEF)
	@echo 
	@echo "* Creating object files:" 

info: 
	@echo
	@echo "target:   $(LIBRARY) $(ROOTMAP)"
	@echo "sources:  $(SOURCES)"
	@echo "headers:  $(HEADERS)"
	@echo "objects:  $(OBJECTS)"
	@echo 
	@echo "compiler: $(CXX)"
	@echo "flags:    $(CXXFLAGS)"
	@echo "libs:     $(LIBS)"
	@echo

clean:
	$(RM) *.o *.d *.d.* *~ *Dict* $(ROOTMAP) $(LIBRARY)

tags:
	ctags --c-kinds=+p $(HEADERS) $(SOURCES)

%.exe:%.C
	$(CXX) $< $(CXXFLAGS) $(LIBS) -L. -l$(LIBNAME) -o $@

.PHONY: all info tags clean

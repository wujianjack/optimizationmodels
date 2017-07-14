# change this if necessary
GRB_VER = 75
MSVC_VER = 2013
PRB_NAME = cutstock
# end change

CXX = cl
CXXFLAGS = /nologo /MD /Ox /EHsc /DNDEBUG 

LD = link
LDFLAGS = /nologo 

INCL = -I$(GUROBI_HOME)\include
ADDINCLFLAGS = 

LIBS = $(GUROBI_HOME)\lib\gurobi$(GRB_VER).lib $(GUROBI_HOME)\lib\gurobi_c++md$(MSVC_VER).lib
ADDLIBS = 

PRB_OBJS = $(PRB_NAME).obj 
PRB_TARGET = $(PRB_NAME).exe 

all: $(PRB_TARGET)

$(PRB_TARGET): $(PRB_OBJS)
	$(LD) $(LDFLAGS) $(LIBS) $(ADDLIBS) $**

.cpp.obj:
	$(CXX) $(CXXFLAGS) -c $(INCL) $(ADDINCLFLAGS) $<
	
clean:
	@del /s /q $(PRB_OBJS) $(PRB_TARGET) > nul
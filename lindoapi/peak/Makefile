# change if necessary
PRB_NAME = peak
LINDOAPI_VERSION = 11
# end change

CC = cl
CFLAGS = /nologo /MT /Ox /EHsc /DNDEBUG /D_LINDO_DLL_

LD = link
LDFALGS = /nologo

INCL = -I"$(LINDOAPI_HOME)include"

LIBS = -libpath:"$(LINDOAPI_HOME)lib\win32" lindo$(LINDOAPI_VERSION)_0.lib

PRB_OBJS = $(PRB_NAME).obj
PRB_TARGET = $(PRB_NAME).exe

all: $(PRB_TARGET)

$(PRB_TARGET): $(PRB_OBJS)
	$(LD) $(LDFALGS) -out:$@ $** $(LIBS)

.c.obj:
	$(CC) $(CFLAGS) $(INCL) -c $<
	
clean:
	@del /s /q $(PRB_OBJS) $(PRB_TARGET) > nul
#Makefile for smwga

CC		= g++
#GDB_FLAGS 	= -Wall -ggdb
CC_OPTS		= -O2
LIBS		= -lm

TARGET		= wave3

SOURCES		= main.cpp random.cpp population.cpp readconf.cpp model.cpp

OBJECTS		= main.o random.o population.o readconf.o model.o

# Implicit rules

#.SUFFIXES:
.SUFFIXES: .cpp

%.o:	%.cpp
	$(CC) -c $(GDB_FLAGS) $(CC_OPTS) $(INCDIR) $<

# Build rules

all:	$(TARGET)

$(TARGET):	$(OBJECTS)
	$(CC) $(GDB_FLAGS) $(CC_OPTS) $(OBJECTS) $(INCDIR) $(LIBDIR) $(LIBS) -o $(TARGET)

depend:
	@makedepend $(SOURCES) 2> /dev/null
	
clean:	
	-rm -f *.o *~ 

# DO NOT DELETE THIS LINE -- make depend depends on it.

main.o: model.h config.h population.h random.h
random.o: random.h
population.o: random.h population.h config.h
readconf.o: model.h config.h population.h random.h
model.o: model.h config.h population.h random.h

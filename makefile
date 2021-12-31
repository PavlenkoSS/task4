CC=g++

LINKER_FLAGS=-c
FLAGS_COMPILE=-O2 -pg -pthread -ffast-math
CPP_FLAGS_COMPATABILITY=-pthread -lm -fsignaling-nans -fstack-protector-all -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format

SOURCES=main.cpp matinit.cpp rotmat.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=mpiths

all: $(OBJECTS)
	$(CC) $(FLAGS_COMPILE) $^ -o $(EXECUTABLE)

.cpp.o:
	$(CC) $(LINKER_FLAGS) $(CPP_FLAGS_COMPATABILITY) $(FLAGS_COMPILE) $< 

clean:
	rm -rf *.o
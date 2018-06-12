CC       = mpicxx.mpich2
CFLAGS   = -O2 -Isrc -Igenericio -Idtk -I/home/dkorytov/proj/hacc/trunk/genericio -g
LDFLAGS  = -Lgenericio/genericio.build/libs -Ldtk/lib -lGenericIO -ldtk
hdf5_opts= /usr/lib/x86_64-linux-gnu/libhdf5_hl_cpp.a /usr/lib/x86_64-linux-gnu/libhdf5_cpp.a /usr/lib/x86_64-linux-gnu/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/libhdf5.a -Bsymbolic-functions -z relro -lpthread -lz -ldl -lstdc++ -lm -lgcc_s -lgcc -lc -lgcc_s -lgcc /usr/lib/gcc/x86_64-linux-gnu/4.8/../../../x86_64-linux-gnu/crtn.o

SOURCES  = $(wildcard src/*.cpp)
OBJECTS  = $(SOURCES:src/%.cpp=obj/%.o)
EXE      = main


${EXE}: ${OBJECTS}
	${CC} ${CFLAGS} -o $@ $^ ${LDFLAGS} ${hdf5_opts}

${OBJECTS}: obj/%.o : src/%.cpp
	${CC} ${CFLAGS} -c -o $@ $<

.PHONY:
clean:
	@rm -f obj/*.o
	@rm -f main


CXX=		g++
CC=			gcc
CXXFLAGS=	-g -O3 -msse4.2 -mpopcnt -fomit-frame-pointer -Wall
CFLAGS=		$(CXXFLAGS)
CPPFLAGS=
INCLUDES=
OBJS=		cov.o
EXE=		cov_cal
LIBS=		-lz -lpthread -lm

ifneq ($(asan),)
	CXXFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.SUFFIXES:.cpp .c .o
.PHONY:all clean depend

.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(EXE)

$(EXE):$(OBJS) main.o
		$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(EXE) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CPPFLAGS) $(DFLAGS) -- *.cpp)

# DO NOT DELETE

cov.o: cov.h ksort.h kvec.h khashl.h kseq.h
main.o: cov.h

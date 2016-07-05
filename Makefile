CC = g++
LD = g++

CFLAGS = -O3
LDFLAGS = -O3

all: spa_run

IPATH = spa/cpp

spa.o: spa.cpp spa.h spa_util.h spa_io.h spa_types.h
	$(CC) -I$(IPATH) -c $(CFLAGS) spa.cpp

spa_util.o: spa_util.cpp spa_util.h
	$(CC) -I$(IPATH) -c $(CFLAGS) spa_util.cpp

spa_io.o: spa_io.cpp spa_types.h
	$(CC) -I$(IPATH) -c $(CFLAGS) spa_io.cpp

spa_run: spa.o spa_util.o spa_io.o
	$(LD) $(LDFLAGS) spa.o spa_util.o spa_io.o -o spa

clean:
	rm -f *.o
	rm spa

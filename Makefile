CC = gcc
LDFLAGS = -lm 
OPT=-O2
CFLAGS = ${OPT} $(INCLUDE) ${INCFLAGS} ${PROF} ${PARALLEL}
test:main.c
#test:metal.c
	$(CC) $(CFLAGS) $? -lgsl -lm -o $@

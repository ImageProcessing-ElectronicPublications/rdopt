PROJECT = rdopt
SRC = src/
CC = gcc
CFLAGS = -O -I$(SRC)
LDFLAGS = -lm -s
SRCS = $(SRC)command.c \
       $(SRC)dct.c \
       $(SRC)histogram.c \
       $(SRC)Image.c \
       $(SRC)jobmethods.c \
       $(SRC)optimize.c \
       $(SRC)utils.c \
       $(SRC)errbpp.c \
       $(SRC)qentry.c \
       $(SRC)lagr_errbpp.c \
       $(SRC)lagr_command.c \
       $(SRC)compress.c \
       $(SRC)usage.c \
       $(SRC)main.c

OBJS = $(SRC)command.o \
       $(SRC)dct.o \
       $(SRC)histogram.o \
       $(SRC)Image.o \
       $(SRC)jobmethods.o \
       $(SRC)optimize.o \
       $(SRC)utils.o \
       $(SRC)errbpp.o \
       $(SRC)qentry.o \
       $(SRC)lagr_errbpp.o \
       $(SRC)lagr_command.o \
       $(SRC)compress.o \
       $(SRC)usage.o \
       $(SRC)main.o

all: $(PROJECT) $(PROJECT)-qttosf

$(PROJECT): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(PROJECT)-qttosf: $(SRC)QTtoSF.c
	$(CC) $(CFLAGS) $^ -o $@ -s

clean:
	rm $(PROJECT) $(PROJECT)-qttosf $(OBJS)

depend:
	makedepend $(SRCS)

$(SRC)command.o: $(SRC)rdopt.h $(SRC)version.h $(SRC)precision.h $(SRC)Image.h $(SRC)opers.h $(SRC)qentry.h
$(SRC)dct.o: $(SRC)rdopt.h $(SRC)version.h $(SRC)precision.h $(SRC)Image.h $(SRC)opers.h $(SRC)CosTable.h
$(SRC)histogram.o: $(SRC)rdopt.h $(SRC)version.h $(SRC)precision.h $(SRC)Image.h $(SRC)opers.h
$(SRC)Image.o: $(SRC)Image.h $(SRC)precision.h
$(SRC)jobmethods.o: $(SRC)rdopt.h $(SRC)version.h $(SRC)precision.h $(SRC)Image.h $(SRC)opers.h
$(SRC)main.o: $(SRC)rdopt.h $(SRC)version.h $(SRC)precision.h $(SRC)Image.h $(SRC)opers.h
$(SRC)optimize.o: $(SRC)rdopt.h $(SRC)version.h $(SRC)precision.h $(SRC)Image.h $(SRC)opers.h
$(SRC)utils.o: $(SRC)rdopt.h $(SRC)version.h $(SRC)precision.h $(SRC)Image.h $(SRC)opers.h
$(SRC)errbpp.o: $(SRC)rdopt.h $(SRC)version.h $(SRC)precision.h $(SRC)Image.h $(SRC)opers.h $(SRC)qentry.h
$(SRC)qentry.o: $(SRC)precision.h $(SRC)qentry.h
$(SRC)lagr_errbpp.o: $(SRC)rdopt.h $(SRC)version.h $(SRC)precision.h $(SRC)Image.h $(SRC)opers.h $(SRC)qentry.h
$(SRC)lagr_command.o: $(SRC)rdopt.h $(SRC)version.h $(SRC)precision.h $(SRC)Image.h $(SRC)opers.h $(SRC)qentry.h
$(SRC)compress.o: $(SRC)rdopt.h $(SRC)version.h $(SRC)precision.h $(SRC)Image.h $(SRC)opers.h

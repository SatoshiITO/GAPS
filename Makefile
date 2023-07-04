CC	= gcc
AR	= ar
CFLAGS	= -g -Wall -O2
LIBS  	= -lhts -lz -lefence
HTSLIB	= $(HOME)/build/htslib2
INCDIR  = -I$(HOME)/.local/include -I$(HTSLIB)
LIBDIR  = -L$(HOME)/.local/lib

PROGRAM	= fastIO

HEADERS = 
OBJS 	= $(SRCS:.c=.o)
SRCS 	= fastIO.c \
	      utils.c



all: $(PROGRAM)

.SUFFIXES: .c .o

.c.o:
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@

$(PROGRAM) : $(OBJS)
	$(CC) $(CFLAGS) $(INCDIR) $(LIBDIR) -o $@ $(OBJS) $(LIBS)

clean:
	rm -rf $(OBJS) *~

realclean: clean
	rm -f $(PROGRAM)

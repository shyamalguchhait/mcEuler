CC = GCC
CFLAGS = -g -Wall -ansi -O3
LD = GCC
LDFLAGS = -lm
RM = /bin/rm -f
OBJS = stok1.o mie.o complex.o nrutil.o array.o
PROG = iquv

# top-level rule, to compile everything.
all: $(PROG)

# rule to link the program
$(PROG): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(PROG)

# rule for file "stok1.o".
stok1R.o: stok1R.c
	$(CC) $(CFLAGS) -c stok1.c

# rule for file "mie.o".
mie.o: mie.c
	$(CC) $(CFLAGS) -c mie.c

# rule for file "nrutil.o".
nrutil.o: nrutil.c
	$(CC) $(CFLAGS) -c nrutil.c

# rule for file "complex.o".
complex.o: complex.c
	$(CC) $(CFLAGS) -c complex.c
	
# rule for file "array.o".
array.o: array.c
	$(CC) $(CFLAGS) -c array.c


# rule for cleaning re-compilable files.
clean:
	$(RM) $(PROG) $(OBJS)

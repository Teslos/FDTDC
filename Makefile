RUNFILE	= quantfdtd	

CC		=	clang

SRC		= 	quantfdtd.c
OBJ		= 	quantfdtd.o

HEADERS	=	

# Optionen fuer Debugger
HOME = /cluster/home/matl/ivast
#CFLAGS 		=	-g -c -I/opt/local/include
# Optionen fuer INDIGO1
# FLAGS		=	-O2 -mips2 -c 
# FLAGS		=	-c -Aa +O4
# FLAGS		=	-inline  -O2 -mips4 -c -WK,-o=5,-so=3,-r=3 -woff all
CFLAGS           = -c -O3 -I/opt/local/include 

LIB 		=  -L/Users/toniivas/Documents/Sources/fftw-3.3/.libs -lc -lm -lfftw3

$(RUNFILE) 	:	$(HEADERS) $(OBJ)
	$(CC) $(OBJ) -o $(RUNFILE) $(LIB) 
#strip $(RUNFILE)

quantfdtd.o : quantfdtd.c
	$(CC) $(CFLAGS) quantfdtd.c
quanttest.o     :	quanttest.c
	$(CC) $(CFLAGS) quanttest.c

clean:
	rm $(OBJ)
	rm $(RUNFILE)


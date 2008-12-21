include ../../Make.config

NAME   = bessel
HERMES = ../../lib/libhermes2d-$(CC)-O-cplx.a
#FLAGS  = -ggdb -DCOMPLEX
FLAGS  = -O3 -DNDEBUG -DCOMPLEX 


$(NAME): main.o bessel.o $(HERMES)
	$(CC) main.o bessel.o -o $@ $(HERMES) $(LIBS) -pthread

main.o: main.cpp ../../Make.config ../../include/*.h
	$(CC) -c $(FLAGS) $< -o $@ -I../../src

bessel.o: bessel.cpp
	$(CC) -c -O3 bessel.cpp

clean:
	rm -f $(NAME) *.o core core.[0-9]*

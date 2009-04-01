include ../../Make.config

NAME   = bessel
HERMES = ../../lib/libhermes2d-$(CC)-O-cplx.a
#FLAGS  = -ggdb -DCOMPLEX
FLAGS  = -O3 -DNDEBUG -DCOMPLEX 


$(NAME): main.o bessel.o $(HERMES)
	$(LD) main.o bessel.o -o $@ $(HERMES) $(LIBS) -pthread 

main.o: main.cpp ../../Make.config
	$(CC) -c $(FLAGS) $< -o $@ -I../../include $(INCLUDES)

bessel.o: bessel.cpp
	$(CC) -c -O3 bessel.cpp

clean:
	rm -f $(NAME) *.o


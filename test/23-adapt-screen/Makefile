include ../../Make.config

NAME   = screen
HERMES = ../../lib/libhermes2d-$(CC)-O-cplx.a
#FLAGS  = -ggdb -DCOMPLEX
FLAGS  = -O3 -DNDEBUG -DCOMPLEX 


$(NAME): main.o fresnl.o $(HERMES)
	g++ main.o fresnl.o -o $@ $(HERMES) $(LIBS) -pthread

main.o: main.cpp ../../Make.config
	g++ -c $(FLAGS) $< -o $@ -I../../src

fresnl.o: fresnl.c ../../Make.config
	gcc -c $< -o $@ 

clean:
	rm -f $(NAME) *.o


include ../../Make.config

HERMES = ../../lib/libhermes2d-$(CC)-O-real.a
FLAGS  = -O3 -DNDEBUG
#FLAGS = -ggdb

all: linview vecview ordview

linview: main.cpp
	$(CC) $(FLAGS) $< -o $@ -DWHATVIEW=ScalarView -DCONTOURS -I../../include $(HERMES) $(LIBS) -pthread

vecview: main.cpp
	$(CC) $(FLAGS) $< -o $@ -DWHATVIEW=VectorView -I../../include $(HERMES) $(LIBS) -pthread

ordview: main.cpp
	$(CC) $(FLAGS) $< -o $@ -DWHATVIEW=OrderView -I../../include $(HERMES) $(LIBS) -pthread


clean:
	rm *.o linview vecview ordview

CXX=mpic++ -Werror -Wall -O3  -ffast-math -mavx -ftree-vectorize -march=native -fopenmp -funroll-loops

CPPFLAGS=-Iinc -g

LDFLAGS=-lboost_program_options -lm

FILES=export forward init memory parse_args shalw
OBJECTS=$(addsuffix .o, $(FILES))
BIN=bin/shalw

all : $(BIN)

$(BIN) : $(addprefix obj/, $(OBJECTS))
	$(CXX) -o $@ $^ $(LDFLAGS)

obj/%.o : src/%.cpp
	$(CXX) -c -o $@ $^ $(CPPFLAGS)

obj/%.o : src/%.c
	$(CXX) -c -o $@ $^ $(CPPFLAGS)

clean :
	rm -f bin/* obj/*

CXX      = /usr/bin/g++
CFLAGS  = -g -DGIT_VERSION=\"$$(git log | head -n1 | cut -f2 -d' ')\" -I /home/mit/csb/fischer/libs/include/ -std=c++11
LDFLAGS = 

OBJ = main.o

jrnf_int: $(OBJ)
	$(CXX) $(CFLAGS) -o jrnf_int $(OBJ) $(LDFLAGS)

clean:
	rm -f $(OBJ); rm -f jrnf_int

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $<




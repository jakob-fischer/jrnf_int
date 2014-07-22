CXX      = /usr/bin/g++
CFLAGS  = -g -DGIT_VERSION=\"$$(git log | head -n1 | cut -f2 -d' ')\" -I /home/mit/csb/fischer/libs/include/ -std=c++11
LDFLAGS = 

OBJ = main.o

prog: $(OBJ)
	$(CXX) $(CFLAGS) -o prog $(OBJ) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $<




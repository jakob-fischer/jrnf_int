CXX      = /usr/bin/g++
CFLAGS  = -g -DGIT_VERSION=\"$$(git log | head -n1 | cut -f2 -d' ')\" -I /home/mit/csb/fischer/libs/include/ -std=c++11
LDFLAGS = 

OBJ = main.o

odeint_rnet: $(OBJ)
	$(CXX) $(CFLAGS) -o odeint_rnet $(OBJ) $(LDFLAGS)

clean:
	rm -f $(OBJ); rm -f odeint_rnet

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $<




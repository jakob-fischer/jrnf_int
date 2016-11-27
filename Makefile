# In case of error: check that g++ is installed (new enough to support c++11) and in the path
# Also boost has to be installed (maybe path to include files has to be given with -I option)
CXX      = g++
CFLAGS  = -g -DGIT_VERSION=\"$$(git log | head -n1 | cut -f2 -d' ')\" -std=c++11 -O4
LDFLAGS = 

OBJ = main.o

jrnf_int: $(OBJ)
	$(CXX) $(CFLAGS) -o jrnf_int $(OBJ) $(LDFLAGS)

clean:
	rm -f $(OBJ); rm -f jrnf_int

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $<




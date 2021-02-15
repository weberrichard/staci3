
-include subdir.mk

CXX=clang++
CXXFLAGS=-g -std=c++11 -Wall -pedantic

# LDFLAGS:   -L/usr/local/lib 
	# -L/usr/lib -L/usr/local/opt -L/usr/local/Cellar/suite-sparse/4.5.3
# CFLAGS:	-I/usr/local/opt/tbb/include 

all: new_staci

new_staci: $(OBJS) $(USER_OBJS)
	$(CXX) $(CXXFLAGS) -I/usr/local/opt/tbb/include -o new_staci $(OBJS) $(USER_OBJS)  -L/usr/local/lib -L/usr/lib -L/usr/local/opt -L/usr/local/Cellar/suite-sparse/4.5.3 -lm -lumfpack #-lga

clean:
	-rm $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) new_staci
	
	# g++ -o new_staci $(OBJS) $(USER_OBJS) -lm -lumfpack -lga

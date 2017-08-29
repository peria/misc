TARGET	= client
SRC	= client.cc connection.cc daihinmin.cc
OBJ	= $(SRC:.cc=.o)

CXX = clang++
LIBS = -lnsl   
CXXFLAGS = -g -O2

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

.cc.o:
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f $(OBJ) $(TARGET)

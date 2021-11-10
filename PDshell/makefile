# The main code and all its relevant files
PdShell: PdShell.cpp PdShell.h user.h

# the compiler: gcc for C program, define as g++ for C++
CC = g++ 

# compiler flags:
#  -Omax  - this flag is used to perform maximum optimization 
#  -g     - this flag adds debugging information to the executable file
#  -Wall  - this flag is used to turn on most compiler warnings
#
#CFLAGS = -g -Wall -Wextra -Wpedantic -Wunused -Wshadow -Wpointer-arith -Wcast-qual -Wmissing-prototypes -Wno-missing-braces -std=c++14 
CFLAGS = -O3 -std=c++14 

TARGET = PdShell

all: ${TARGET}

$(TARGET): $(TARGET).cpp
					 $(CC) $(CFLAGS) -o $(TARGET) $(TARGET).cpp

clean: 
	$(RM) $(TARGET)

CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall
TARGET = ClassicalXX
SRC = main.cpp
CSV = evolution.csv

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET) $(CSV)

run: clean
	make && ./simu && python3 plot.py


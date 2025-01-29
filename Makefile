TARGET = prog
LIBS = -lm -lc
CC = g++
CFLAGS = -flto -g -std=c++23 -Wall -O3
LFLAGS = -flto -g -O3

# Source directories
SRC_DIRS = .

# Create obj directory path for each source directory
OBJ_DIR = obj
OBJ_DIRS = $(foreach dir,$(SRC_DIRS),$(OBJ_DIR)/$(dir))

SOURCES = $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.cpp))

OBJECTS = $(SOURCES:%.cpp=$(OBJ_DIR)/%.o)

INCLUDES = $(foreach dir,$(SRC_DIRS),-I$(dir))

.PHONY: default all clean

default: $(TARGET)
all: default

# Create necessary directories before compilation
$(OBJ_DIR):
	mkdir -p $(OBJ_DIRS)

$(OBJ_DIR)/%.o: %.cpp | $(OBJ_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -Wall $(LIBS) -o $@

clean:
	-rm -rf $(OBJ_DIR)
	-rm -f $(TARGET)

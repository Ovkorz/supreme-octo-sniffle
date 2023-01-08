CC = g++
LD = g++

VERBOSE		= false
CFLAG       = -Wall
PROG_NAME   = prog

ifeq ($(strip $(VERBOSE)), true)
	CFLAG	+= -DVERBOSE
endif

INC = -I ./include
SRC_DIR     = ./source
BUILD_DIR   = ./build
BIN_DIR     = ./bin

SRC_LIST = $(wildcard $(SRC_DIR)/*.cpp)
OBJ_LIST = $(addprefix $(BUILD_DIR)/, $(notdir $(SRC_LIST:.cpp=.o)))

.PHONY: all clean $(PROG_NAME) compile

all: $(PROG_NAME)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) -c $^ $(CFLAG) -o $@ $(INC)

$(PROG_NAME): $(OBJ_LIST)
	$(LD) $(OBJ_LIST) -o $(BIN_DIR)/$@

clean:
	rm -f $(BIN_DIR)/$(PROG_NAME) $(BUILD_DIR)/*.o

# kompilator c
CCOMP = g++

# konsolidator
LOADER = g++

# opcje optymalizacji:
# wersja do debugowania
OPT = -g -DDEBUG
# wersja zoptymalizowana do mierzenia czasu
# OPT = -O3

# pliki naglowkowe
INC = -I ./include
# ALG_INC = -I ./include

# biblioteki
# LIB = -L../pomiar_czasu -lpomiar_czasu -lm

# ./bin/optimize : bin/main.o bin/opt_alg.o	
# 	$(LOADER) $(OPT) ./bin/main.o .bin/opt_alg.o -o .bin/optimize

# ./bin/main.o : main.cpp
# 	$(CCOMP) -c $(OPT) main.cpp -o ./bin/main.o $(INC)

# ./bin/opt_alg.o : ./alg/opt_alg.cpp ./alg/include/opt_alg.h
# 	$(CCOMP) -c $(OPT) ./alg/opt_alg.cpp -o ./bin/opt_alg.o $(ALG_INC)

CFLAG       = -Wall
PROG_NAME   = prog

SRC_DIR     = ./source
BUILD_DIR   = ./build
BIN_DIR     = ./bin
SRC_LIST = $(wildcard $(SRC_DIR)/*.cpp)
OBJ_LIST = $(BUILD_DIR)/$(notdir $(SRC_LIST:.cpp=.o))

.PHONY: all clean $(PROG_NAME) compile

all: $(PROG_NAME)

compile: 
	$(CC) $(SRC_LIST) -o $(BUILD_DIR)/$(OBJ_LIST) $(INC)

$(PROG_NAME): compile
	$(LD) $(OBJ_LIST) -o $(BIN_DIR)/$@

clean:
	rm -f $(BIN_DIR)/$(PROG_NAME) $(BUILD_DIR)/*.o

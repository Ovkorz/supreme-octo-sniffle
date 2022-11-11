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
INC = -I ./alg/include
ALG_INC = -I ./include

# biblioteki
# LIB = -L../pomiar_czasu -lpomiar_czasu -lm

./bin/optimize : bin/main.o bin/opt_alg.o	
	$(LOADER) $(OPT) ./bin/main.o .bin/opt_alg.o -o .bin/optimize

./bin/main.o : main.cpp
	$(CCOMP) -c $(OPT) main.cpp -o ./bin/main.o $(INC)

./bin/opt_alg.o : ./alg/opt_alg.cpp ./alg/include/opt_alg.h
	$(CCOMP) -c $(OPT) ./alg/opt_alg.cpp -o ./bin/opt_alg.o $(ALG_INC)
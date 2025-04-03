# Nom de l'exécutable
EXE = main

# Répertoires
BIN = ./
OBJ = obj
SRC = src
INC = include
IMGUI = imgui
BACKENDS = imgui/backends

# Compilateur
CC = g++
LD = g++

# Flags de compilation
CXXFLAGS = -Wall -O0 -g -I $(INC) -I $(IMGUI) -I $(BACKENDS) 
DEPFLAGS = -MMD -MP

# Flags de l'éditeur de liens
LDFLAGS = -lglfw -lGL -ldl -lpthread -lX11

# Détection des fichiers sources
SOURCES := $(wildcard $(SRC)/*.cpp) $(wildcard $(IMGUI)/*.cpp) $(wildcard $(BACKENDS)/*.cpp)

# Génération des fichiers objets correspondants
OBJECTS := $(patsubst %.cpp, $(OBJ)/%.o, $(SOURCES))

# Inclusion des fichiers de dépendances générés
DEPENDS := $(OBJECTS:.o=.d)

# Compilation des sources
$(OBJ)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CC) $(DEPFLAGS) $(CXXFLAGS) -c $< -o $@

# Règle principale
.PHONY: all
all: $(BIN)/$(EXE)

$(BIN)/$(EXE): $(OBJECTS)
	$(LD) $(OBJECTS) $(LDFLAGS) -o $@

# Exécution
.PHONY: run
run: $(BIN)/$(EXE)
	./$(BIN)/$(EXE)

# Nettoyage
.PHONY: clean
clean:
	$(RM) $(OBJECTS) $(DEPENDS) $(BIN)/$(EXE)

# Nettoyage complet
.PHONY: reset
reset:
	$(RM) -r $(OBJ) $(BIN)

# Inclusion des dépendances
-include $(DEPENDS)

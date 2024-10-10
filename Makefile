BIN   := ./bin
INC   := ./include
SRC   := ./src
SRC_C := ./src-corrections

CXX    := g++ -std=c++11
CFLAGS := -Wall -g -O3

ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTINCDIR  := $(shell root-config --incdir)
ROOTLIBS    := $(shell root-config --libs) -lEG

all: ${BIN}/create_e2c_ntuple ${BIN}/create_e2c_dtrmatchntuple

${BIN}/create_e2c_ntuple: ${SRC}/create_e2c_ntuple.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_e2c_ntuple.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_e2c_ntuple

${BIN}/create_e2c_dtrmatchntuple: ${SRC_C}/create_e2c_dtrmatchntuple.cpp
	${CXX} ${ROOTCFLAGS} ${SRC_C}/create_e2c_dtrmatchntuple.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_e2c_dtrmatchntuple

clean:
	rm ${BIN}/*

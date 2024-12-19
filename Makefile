BIN   := ./bin
INC   := ./include
SRC   := ./src
SRC_C := ./src-analysis

CXX    := g++ -std=c++11
CFLAGS := -Wall -g -O3

ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTINCDIR  := $(shell root-config --incdir)
ROOTLIBS    := $(shell root-config --libs) -lEG

all: ${BIN}/create_e2c_ntuple ${BIN}/create_e2c_purityntuple ${BIN}/create_e2c_corrntuple ${BIN}/create_e2c_efficiencyntuple 

${BIN}/create_e2c_ntuple: ${SRC}/create_e2c_ntuple.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_e2c_ntuple.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_e2c_ntuple

${BIN}/create_e2c_purityntuple: ${SRC_C}/create_e2c_purityntuple.cpp
	${CXX} ${ROOTCFLAGS} ${SRC_C}/create_e2c_purityntuple.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_e2c_purityntuple

${BIN}/create_e2c_corrntuple: ${SRC_C}/create_e2c_corrntuple.cpp
	${CXX} ${ROOTCFLAGS} ${SRC_C}/create_e2c_corrntuple.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_e2c_corrntuple

${BIN}/create_e2c_efficiencyntuple: ${SRC_C}/create_e2c_efficiencyntuple.cpp
	${CXX} ${ROOTCFLAGS} ${SRC_C}/create_e2c_efficiencyntuple.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_e2c_efficiencyntuple

clean:
	rm ${BIN}/*

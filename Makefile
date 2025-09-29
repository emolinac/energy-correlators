BIN := ./bin
INC := ./include
SRC := ./src

CXX    := g++ -std=c++11
CFLAGS := -Wall -g -O3

ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTINCDIR  := $(shell root-config --incdir)
ROOTLIBS    := $(shell root-config --libs) -lEG

all: ${BIN}/create_jet_purityntuple ${BIN}/create_jet_efficiencyntuple \
     ${BIN}/create_jet_purityntuple_ct ${BIN}/create_jet_efficiencyntuple_ct \
     ${BIN}/create_correec_histo3dpaircorr_rl_jetpt_weightpt \
     ${BIN}/create_correec_histo3dpaircorr_rl_jetpt_weightpt_ct \
     ${BIN}/create_eec_paircorrectionsntuple ${BIN}/create_eec_paircorrectionsntuple_ct 

${BIN}/create_jet_purityntuple: ${SRC}/create_jet_purityntuple.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_jet_purityntuple.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_jet_purityntuple

${BIN}/create_jet_efficiencyntuple: ${SRC}/create_jet_efficiencyntuple.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_jet_efficiencyntuple.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_jet_efficiencyntuple

${BIN}/create_jet_purityntuple_ct: ${SRC}/create_jet_purityntuple_ct.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_jet_purityntuple_ct.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_jet_purityntuple_ct

${BIN}/create_jet_efficiencyntuple_ct: ${SRC}/create_jet_efficiencyntuple_ct.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_jet_efficiencyntuple_ct.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_jet_efficiencyntuple_ct

${BIN}/create_correec_histo3dpaircorr_rl_jetpt_weightpt: ${SRC}/create_correec_histo3dpaircorr_rl_jetpt_weightpt.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_correec_histo3dpaircorr_rl_jetpt_weightpt.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_correec_histo3dpaircorr_rl_jetpt_weightpt

${BIN}/create_correec_histo3dpaircorr_rl_jetpt_weightpt_ct: ${SRC}/create_correec_histo3dpaircorr_rl_jetpt_weightpt_ct.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_correec_histo3dpaircorr_rl_jetpt_weightpt_ct.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_correec_histo3dpaircorr_rl_jetpt_weightpt_ct

${BIN}/create_eec_paircorrectionsntuple: ${SRC}/create_eec_paircorrectionsntuple.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_eec_paircorrectionsntuple.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_eec_paircorrectionsntuple

${BIN}/create_eec_paircorrectionsntuple_ct: ${SRC}/create_eec_paircorrectionsntuple_ct.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_eec_paircorrectionsntuple_ct.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_eec_paircorrectionsntuple_ct

${BIN}/create_eec_corrntuple: ${SRC}/create_eec_corrntuple.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_eec_corrntuple.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_eec_corrntuple

clean:
	rm ${BIN}/*

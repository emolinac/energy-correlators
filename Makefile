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
	 ${BIN}/create_eec_corrntuple_paircorr \
	 ${BIN}/create_eec_corrntuple_paircorr_3d \
	 ${BIN}/create_eec_corrntuple_paircorr_ct \
	 ${BIN}/create_eec_corrntuple_paircorr_jes ${BIN}/create_eec_corrntuple_paircorr_jer \
	 ${BIN}/create_eec_mc_ntuple ${BIN}/create_hadron_ntuple ${BIN}/create_jes_jer_ntuple \
	 ${BIN}/create_eec_paircorrectionsntuple ${BIN}/create_eec_paircorrectionsntuple_ct \
	 ${BIN}/create_eec_corrntuple ${BIN}/create_eec_hadroncorrectionsntuple

${BIN}/create_jet_purityntuple: ${SRC}/create_jet_purityntuple.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_jet_purityntuple.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_jet_purityntuple

${BIN}/create_jet_efficiencyntuple: ${SRC}/create_jet_efficiencyntuple.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_jet_efficiencyntuple.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_jet_efficiencyntuple

${BIN}/create_eec_paircorrectionsntuple: ${SRC}/create_eec_paircorrectionsntuple.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_eec_paircorrectionsntuple.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_eec_paircorrectionsntuple

${BIN}/create_eec_paircorrectionsntuple_ct: ${SRC}/create_eec_paircorrectionsntuple_ct.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_eec_paircorrectionsntuple_ct.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_eec_paircorrectionsntuple_ct

${BIN}/create_eec_corrntuple: ${SRC}/create_eec_corrntuple.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_eec_corrntuple.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_eec_corrntuple

${BIN}/create_eec_hadroncorrectionsntuple: ${SRC}/create_eec_hadroncorrectionsntuple.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_eec_hadroncorrectionsntuple.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_eec_hadroncorrectionsntuple

${BIN}/create_eec_corrntuple_paircorr: ${SRC}/create_eec_corrntuple_paircorr.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_eec_corrntuple_paircorr.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_eec_corrntuple_paircorr

${BIN}/create_eec_corrntuple_paircorr_3d: ${SRC}/create_eec_corrntuple_paircorr_3d.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_eec_corrntuple_paircorr_3d.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_eec_corrntuple_paircorr_3d

${BIN}/create_eec_corrntuple_paircorr_ct: ${SRC}/create_eec_corrntuple_paircorr_ct.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_eec_corrntuple_paircorr_ct.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_eec_corrntuple_paircorr_ct

${BIN}/create_eec_corrntuple_paircorr_jes: ${SRC}/create_eec_corrntuple_paircorr_jes.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_eec_corrntuple_paircorr_jes.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_eec_corrntuple_paircorr_jes

${BIN}/create_eec_corrntuple_paircorr_jer: ${SRC}/create_eec_corrntuple_paircorr_jer.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_eec_corrntuple_paircorr_jer.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_eec_corrntuple_paircorr_jer

${BIN}/create_eec_mc_ntuple: ${SRC}/create_eec_mc_ntuple.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_eec_mc_ntuple.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_eec_mc_ntuple

${BIN}/create_eec_mc_at_ntuple: ${SRC}/create_eec_mc_at_ntuple.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_eec_mc_at_ntuple.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_eec_mc_at_ntuple

${BIN}/create_hadron_ntuple: ${SRC}/create_hadron_ntuple.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_hadron_ntuple.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_hadron_ntuple

${BIN}/create_jes_jer_ntuple: ${SRC}/create_jes_jer_ntuple.cpp
	${CXX} ${ROOTCFLAGS} ${SRC}/create_jes_jer_ntuple.cpp -I${INC} ${ROOTLIBS} -o ${BIN}/create_jes_jer_ntuple

clean:
	rm ${BIN}/*

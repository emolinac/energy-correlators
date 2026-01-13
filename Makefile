BIN := ./bin
INC := ./include
SRC := ./src

CXX    := g++ -std=c++11 -march=native
CFLAGS := -Wall -g -O3 

ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTINCDIR  := $(shell root-config --incdir)
ROOTLIBS    := $(shell root-config --libs --glibs) -lEG

CUSTOM_FUNCTIONS := $(wildcard ${INC}/*.cpp)

all: ${BIN}/create_jet_ntuple_truth2reco_match \
     ${BIN}/create_jet_ntuple_truth2reco_match_ct \
     ${BIN}/create_eec_mc_ntuple \
     ${BIN}/create_eec_28r1_mc_ntuple \
     ${BIN}/create_eec_28r2_mc_ntuple \
     ${BIN}/create_jes_jer_ntuple \
     ${BIN}/create_hadron_jet_ntuples \
     ${BIN}/create_correec_histo3dpaircorr_rl_jetpt_weightpt \
     ${BIN}/create_correec_histo3dpaircorr_rl_jetpt_weightpt_ct \
     ${BIN}/create_ntuple_reco2truth_match ${BIN}/create_ntuple_reco2truth_match_ct \
     ${BIN}/create_ntuple_reco2truth_singlehadron_match

${BIN}/create_eec_mc_ntuple: ${SRC}/create_eec_mc_ntuple.cpp
	${CXX} ${ROOTCFLAGS} -I${INC} ${SRC}/create_eec_mc_ntuple.cpp ${CUSTOM_FUNCTIONS} ${ROOTLIBS} -o ${BIN}/create_eec_mc_ntuple

${BIN}/create_eec_28r1_mc_ntuple: ${SRC}/create_eec_28r1_mc_ntuple.cpp
	${CXX} ${ROOTCFLAGS} -I${INC} ${SRC}/create_eec_28r1_mc_ntuple.cpp ${CUSTOM_FUNCTIONS} ${ROOTLIBS} -o ${BIN}/create_eec_28r1_mc_ntuple

${BIN}/create_eec_28r2_mc_ntuple: ${SRC}/create_eec_28r2_mc_ntuple.cpp
	${CXX} ${ROOTCFLAGS} -I${INC} ${SRC}/create_eec_28r2_mc_ntuple.cpp ${CUSTOM_FUNCTIONS} ${ROOTLIBS} -o ${BIN}/create_eec_28r2_mc_ntuple

${BIN}/create_hadron_jet_ntuples: ${SRC}/create_hadron_jet_ntuples.cpp
	${CXX} ${ROOTCFLAGS} -I${INC} ${SRC}/create_hadron_jet_ntuples.cpp ${CUSTOM_FUNCTIONS} ${ROOTLIBS} -o ${BIN}/create_hadron_jet_ntuples

${BIN}/create_jes_jer_ntuple: ${SRC}/create_jes_jer_ntuple.cpp
	${CXX} ${ROOTCFLAGS} -I${INC} ${SRC}/create_jes_jer_ntuple.cpp ${CUSTOM_FUNCTIONS} ${ROOTLIBS} -o ${BIN}/create_jes_jer_ntuple

${BIN}/create_jet_ntuple_truth2reco_match: ${SRC}/create_jet_ntuple_truth2reco_match.cpp
	${CXX} ${ROOTCFLAGS} -I${INC} ${SRC}/create_jet_ntuple_truth2reco_match.cpp ${CUSTOM_FUNCTIONS} ${ROOTLIBS}  -o ${BIN}/create_jet_ntuple_truth2reco_match

${BIN}/create_jet_ntuple_truth2reco_match_ct: ${SRC}/create_jet_ntuple_truth2reco_match_ct.cpp
	${CXX} ${ROOTCFLAGS} -I${INC} ${SRC}/create_jet_ntuple_truth2reco_match_ct.cpp ${CUSTOM_FUNCTIONS} ${ROOTLIBS} -o ${BIN}/create_jet_ntuple_truth2reco_match_ct

${BIN}/create_correec_histo3dpaircorr_rl_jetpt_weightpt: ${SRC}/create_correec_histo3dpaircorr_rl_jetpt_weightpt.cpp
	${CXX} ${ROOTCFLAGS} -I${INC} ${SRC}/create_correec_histo3dpaircorr_rl_jetpt_weightpt.cpp ${CUSTOM_FUNCTIONS} ${ROOTLIBS} -o ${BIN}/create_correec_histo3dpaircorr_rl_jetpt_weightpt

${BIN}/create_correec_histo3dpaircorr_rl_jetpt_weightpt_ct: ${SRC}/create_correec_histo3dpaircorr_rl_jetpt_weightpt_ct.cpp
	${CXX} ${ROOTCFLAGS} -I${INC} ${SRC}/create_correec_histo3dpaircorr_rl_jetpt_weightpt_ct.cpp ${CUSTOM_FUNCTIONS} ${ROOTLIBS} -o ${BIN}/create_correec_histo3dpaircorr_rl_jetpt_weightpt_ct

${BIN}/create_ntuple_reco2truth_match: ${SRC}/create_ntuple_reco2truth_match.cpp
	${CXX} ${ROOTCFLAGS} -I${INC} ${SRC}/create_ntuple_reco2truth_match.cpp ${CUSTOM_FUNCTIONS} ${ROOTLIBS} -o ${BIN}/create_ntuple_reco2truth_match

${BIN}/create_ntuple_reco2truth_singlehadron_match: ${SRC}/create_ntuple_reco2truth_singlehadron_match.cpp
	${CXX} ${ROOTCFLAGS} -I${INC} ${SRC}/create_ntuple_reco2truth_singlehadron_match.cpp ${CUSTOM_FUNCTIONS} ${ROOTLIBS} -o ${BIN}/create_ntuple_reco2truth_singlehadron_match

${BIN}/create_ntuple_reco2truth_match_ct: ${SRC}/create_ntuple_reco2truth_match_ct.cpp
	${CXX} ${ROOTCFLAGS} -I${INC} ${SRC}/create_ntuple_reco2truth_match_ct.cpp ${CUSTOM_FUNCTIONS} ${ROOTLIBS} -o ${BIN}/create_ntuple_reco2truth_match_ct

${BIN}/create_eec_corrntuple: ${SRC}/create_eec_corrntuple.cpp
	${CXX} ${ROOTCFLAGS} -I${INC} ${SRC}/create_eec_corrntuple.cpp ${CUSTOM_FUNCTIONS} ${ROOTLIBS} -o ${BIN}/create_eec_corrntuple

clean:
	rm ${BIN}/*

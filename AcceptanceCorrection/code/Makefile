CXX = g++
CXXFLAGS = -O3 -Wall -fPIC -std=c++11 -fpermissive


ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --evelibs)
CXXFLAGS += $(ROOTFLAGS)# -mavx
LDLIBS   += $(ROOTLIBS) -lEGPythia6

ana: McEtaF_QA19.C
	$(CXX) $(CXXFLAGS) -o ana McEtaF_QA19.C $(LDLIBS)

clean:
	rm ana 


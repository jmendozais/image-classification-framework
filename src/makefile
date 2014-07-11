CC		  := g++
LD      := g++
CFLAGS  := -lboost_system-mt -lboost_filesystem-mt -fopenmp `pkg-config --cflags --libs opencv` -MP -MD -c  
LDFLAGS := -lboost_system-mt -lboost_filesystem-mt -fopenmp `pkg-config --cflags --libs opencv` -MP -MD 
SOURCES	  := $(shell echo new/*.cpp classification/*.cpp feature_extraction/*.cpp feature_selection/*.cpp ipcore/*.cpp)
HEADERS	  := $(shell echo new/*.h ./*.h classification/*.h feature_extraction/*.h feature_selection/*.h image_transformation/*.h ipcore/*.h)
OBJECTS	  := $(SOURCES:.cpp=.o)
TESTS   := $(shell echo ./*.cpp)

all: $(SOURCES) $(TESTS)

$(TESTS) : $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $@ -o $@.o

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(OBJECTS)

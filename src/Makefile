include Makefile_common

TARGET=main
PROG_NAME=vc_deg
SUBDIRS=util solv
LIB_CPPS=$(shell ls -f $(addsuffix /*.cpp,$(SUBDIRS)))
LIB_OS=$(LIB_CPPS:.cpp=.o)

all: $(SUBDIRS)
	g++ $(CFLAGS) -std=c++0x -Wall -static  ${LIB_OS} ${TARGET}.cpp -o ${PROG_NAME} 


$(SUBDIRS):
	make -s -C $@

tests:
	cd tests && { ./testing ; cd .. ; }

clean:
	rm -f $(shell find -name "*.o") ${PROG_NAME}

.PHONY: $(SUBDIRS) tests clean

ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

export LIBCUCKOO = -I${ROOT_DIR}/lib/libcuckoo
export INTERVAL_TREE = -I${ROOT_DIR}/lib/interval_tree
export LEMON = -I${ROOT_DIR}/lib/lemon
export BIN_DIR = ${ROOT_DIR}/bin
export MINIMAP2_DIR = ${ROOT_DIR}/lib/minimap2

export CXXFLAGS += ${LIBCUCKOO} ${INTERVAL_TREE} ${LEMON} -I${MINIMAP2_DIR} -I ${ROOT_DIR}
export LDFLAGS += -lz -L${MINIMAP2_DIR} -lminimap2

.PHONY: clean all profile debug minimap2

.DEFAULT_GOAL := all


${BIN_DIR}/flye-minimap2:
	make -C ${MINIMAP2_DIR}
	cp ${MINIMAP2_DIR}/minimap2 ${BIN_DIR}/flye-minimap2

minimap2: ${BIN_DIR}/flye-minimap2

all: minimap2
	make release -C src
profile: minimap2
	make profile -C src
debug: minimap2
	make debug -C src
clean:
	make clean -C src
	make clean -C ${MINIMAP2_DIR}
	rm ${BIN_DIR}/flye-minimap2

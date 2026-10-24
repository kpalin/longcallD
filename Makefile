GCC_CHECK := $(shell gcc --version | head -n 1 | grep -i "clang")

# Check if the OS is macOS or linux
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

ifeq ($(GCC_CHECK),) # gcc
	CXXFLAGS = -std=c++11
else # clang
	CXXFLAGS = -std=c++11 -stdlib=libc++
endif

# add -fno-tree-vectorize to avoid certain vectorization errors in O3 optimization
# right now, we are using -O3 for the best performance, and no vectorization errors were found
EXTRA_FLAGS = -Wall -Wno-misleading-indentation -Wno-unused-function #-Wno-unused-variable -Wno-alloc-size-larger-than

# Define the version number
VERSION=0.0.11
# Get the Git commit hash
GIT_COMMIT := $(shell git rev-parse --short HEAD 2> /dev/null)
ifneq ($(GIT_COMMIT),)
	LONGCALLD_VERSION = $(VERSION)-$(GIT_COMMIT)
else
	LONGCALLD_VERSION = $(VERSION)
endif

HTSLIB_DIR  = ./htslib
HTSLIB      = $(HTSLIB_DIR)/libhts.a

EDLIB_DIR     = ./edlib
EDLIB_INC_DIR = ./edlib/include
EDLIB         = $(EDLIB_DIR)/src/edlib.o

ABPOA_DIR   = ./abPOA
ABPOA_LIB   = ./abPOA/lib/libabpoa.a
ABPOA_INC_DIR = $(ABPOA_DIR)/include

WFA2_DIR    = ./WFA2-lib
WFA2_LIB    = $(WFA2_DIR)/lib/libwfa.a

LIB_PATH    =
LIB         = $(HTSLIB) $(ABPOA_LIB) $(WFA2_LIB) $(LIB_PATH) -lm -lz -lpthread -llzma -lbz2 -lcurl
INCLUDE     = -I $(HTSLIB_DIR) -I $(EDLIB_INC_DIR) -I $(ABPOA_INC_DIR) -I $(WFA2_DIR)

# Try linking against libdeflate
ifeq ($(shell echo "int main() {return 0;}" | ${CC} -x c - $(LIB_PATH) -ldeflate >/dev/null 2>&1 && echo "yes"),yes)
	LIB += -ldeflate
endif

ifeq ($(UNAME_S),Linux) # Linux
	LIB += -lcrypto
	ifneq ($(portable),)
		LIB += -static-libgcc -static-libstdc++
		ifneq ($(opt_lib),)
			LIB = $(HTSLIB) $(ABPOA_LIB) $(WFA2_LIB) -static-libgcc -static-libstdc++ -L${opt_lib} -lm -lz -lpthread -llzma -lbz2 -lcurl -lssl -lcrypto -lssh2 -ldl -Wl,-Bstatic -ldeflate -Wl,-Bdynamic -lzstd
		else
			ifneq ($(OPT_LIB),)
				LIB = $(HTSLIB) $(ABPOA_LIB) $(WFA2_LIB) -static-libgcc -static-libstdc++ -L${OPT_LIB} -lm -lz -lpthread -llzma -lbz2 -lcurl -lssl -lcrypto -lssh2 -ldl -Wl,-Bstatic -ldeflate -Wl,-Bdynamic -lzstd
			endif
		endif
	endif
endif

# for debug
ifneq ($(debug),)
	EXTRA_FLAGS  += -D __DEBUG__
endif

ABPOA_GDB_LIB = ./abPOA/lib/libabpoa_sse41.a
ABPOA_NOR_LIB = ./abPOA/lib/libabpoa.a
WFA_GDB_LIB   = ./WFA2-lib/lib/libwfa_gdb.a
WFA_NOR_LIB   = ./WFA2-lib/lib/libwfa.a

# for gdb
ifneq ($(gdb),)
	OPT_FLAGS = -O0 -g
	ABPOA_LIB = $(ABPOA_GDB_LIB)
	WFA2_LIB  = $(WFA_GDB_LIB)
else
	OPT_FLAGS = -O3
	ABPOA_LIB = $(ABPOA_NOR_LIB)
	WFA2_LIB  = $(WFA_NOR_LIB)
endif

# for sanitizers
#   make san=1       AddressSanitizer + UndefinedBehaviorSanitizer
#   make san=thread  ThreadSanitizer
#
# Only longcallD's own objects are instrumented here. htslib, abPOA and WFA2
# link in as prebuilt static libraries and stay uninstrumented; ASan's
# allocator interceptors still catch heap overflow and use-after-free on any
# allocation regardless of which object made it. To instrument htslib as well:
#   make clean_hts && make san=1 HTSLIB_CFLAGS="-g -O1 -fsanitize=address"
#
# We deliberately keep the non-gdb abPOA/WFA2 libraries: WFA2-lib provides no
# libwfa_gdb.a target, so $(WFA_GDB_LIB) is a dangling path and selecting it
# would break the link.
ifeq ($(san),thread)
	OPT_FLAGS = -O1 -g -fno-omit-frame-pointer
	SAN_FLAGS = -fsanitize=thread
else ifneq ($(san),)
	OPT_FLAGS = -O1 -g -fno-omit-frame-pointer
	SAN_FLAGS = -fsanitize=address,undefined -fno-sanitize-recover=all
endif

# -MMD -MP makes the compiler emit a .d file of real header dependencies
# next to each .o, which -include pulls in below.
CFLAGS = $(OPT_FLAGS) $(EXTRA_FLAGS) $(SAN_FLAGS) -MMD -MP -DLONGCALLD_VERSION=\"$(LONGCALLD_VERSION)\"

BIN_LDFLAGS =

ifeq ($(UNAME_S),Linux)
    # ThreadSanitizer requires a position-independent executable: with -no-pie
    # its shadow memory mapping collides with the executable's own and it
    # aborts at startup. ASan and UBSan work either way.
    ifneq ($(san),thread)
        BIN_LDFLAGS += -no-pie
    endif
endif

BIN_LDFLAGS += $(SAN_FLAGS)


HTSLIB_CFLAGS = -g -O2
ifeq ($(UNAME_S),Linux)
    ifeq ($(UNAME_M),x86_64)
        HTSLIB_CFLAGS += -mssse3
    endif
endif

# for gprof
ifneq ($(pg),)
	OPT_FLAGS = -O3
	PG_FLAG   = -pg
	CFLAGS   += -pg
endif

ifneq ($(PREFIX),)
	OUT_PRE_DIR = $(PREFIX)
else
	OUT_PRE_DIR = .
endif

BIN_DIR = $(OUT_PRE_DIR)/bin
INC_DIR = ./include
SRC_DIR = ./src

HTS_ALL = hts_all
EDLIB_ALL = edlib_all
ABPOA_ALL = abpoa_all
WFA2_ALL = wfa2_all
CSOURCE    = $(wildcard ${SRC_DIR}/*.c)
CPPSOURCE  = $(wildcard $(SRC_DIR)/*.cpp)
CPPSOURCE += $(EDLIB_DIR)/src/edlib.cpp

OBJS    = $(CSOURCE:.c=.o) $(CPPSOURCE:.cpp=.o)
BIN     = $(BIN_DIR)/longcallD
ifneq ($(gdb),)
	BIN = $(BIN_DIR)/gdb_longcallD
endif

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

# `all` must stay the default goal: the FORCE target below would otherwise
# become it, since make picks the first target that does not begin with a dot.
.DEFAULT_GOAL := all

# Switching between a normal build, san=1 and san=thread changes which runtime
# every object needs, but make only rebuilds what the *sources* touched -- so a
# plain `make` after a sanitizer run fails at link time with undefined __tsan_*
# references, and a partially-rebuilt tree can segfault on an ABI mismatch.
#
# Keep the current mode in a stamp file, rewritten at parse time and only when
# it actually changes, so its mtime moves exactly then. Every object depends on
# it, so make invalidates them by its normal rules. (Deleting objects from a
# recipe instead does not work: make has already stat'ed them and will link the
# stale set.)
SAN_MODE := $(if $(san),$(san),none)
$(shell [ "$$(cat .build-mode 2>/dev/null)" = "$(SAN_MODE)" ] || printf '%s' "$(SAN_MODE)" > .build-mode)

$(OBJS): .build-mode

# The commit hash reaches the binary as a -D on CFLAGS, which the compiler's
# dependency tracking cannot see. So an incremental build after a commit left
# whichever objects happened not to need recompiling stamped with the previous
# commit, and bin/longcallD reported a version it was not built from -- in a
# VCF's ##source line and a BAM's @PG VN, the one field that records which
# build produced a result. Same stamp-file idiom as .build-mode above, and the
# object list comes from grep so a new user of the macro is covered by itself.
$(shell [ "$$(cat .build-version 2>/dev/null)" = "$(LONGCALLD_VERSION)" ] || printf '%s' "$(LONGCALLD_VERSION)" > .build-version)
VERSION_OBJS := $(patsubst %.c,%.o,$(shell grep -l LONGCALLD_VERSION $(SRC_DIR)/*.c 2>/dev/null))
$(VERSION_OBJS): .build-version

all: $(HTS_ALL) $(EDLIB) $(ABPOA_LIB) $(WFA2_LIB) $(BIN)

$(HTS_ALL): $(HTSLIB)

$(HTSLIB): $(HTSLIB_DIR)/configure.ac
	cd $(HTSLIB_DIR); autoreconf -i; ./configure CC=${CC} CFLAGS="$(HTSLIB_CFLAGS)"; make CC=${CC} CFLAGS="$(HTSLIB_CFLAGS)"

$(EDLIB): $(EDLIB_DIR)/src/edlib.cpp $(EDLIB_DIR)/include/edlib.h
	$(CXX) $(CFLAGS) $(CXXFLAGS) -c $< $(INCLUDE) -o $@
$(EDLIB_ALL): $(EDLIB)

$(ABPOA_GDB_LIB): 
	cd $(ABPOA_DIR); make libabpoa gdb=1 sse41=1
$(ABPOA_NOR_LIB):
	cd $(ABPOA_DIR); make libabpoa

$(ABPOA_ALL): $(ABPOA_LIB)

$(WFA2_LIB):
	cd $(WFA2_DIR); make setup lib_wfa CC_FLAGS="$(CC_FLAGS) -O3" 
$(WFA2_ALL): $(WFA2_LIB)


$(BIN): $(OBJS) $(ABPOA_LIB) $(HTSLIB) $(WFA2_LIB)
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CXX) $(BIN_LDFLAGS) $(OBJS) -o $@ $(LIB) $(PG_FLAG)
# 	$(CC) $(OBJS) -o $@ $(LIB) $(PG_FLAG)

$(SRC_DIR)/align.o: $(SRC_DIR)/align.c
	$(CC) -c -DUSE_SIMDE -DSIMDE_ENABLE_NATIVE_ALIASES $(CFLAGS) $< $(INCLUDE) -o $@

# Header dependencies come from the compiler (-MMD -MP above), not from a
# hand-maintained list. The list that used to live here was incomplete and
# partly wrong: assign_hap.o had no rule at all, vcf_utils.o did not mention
# bam_utils.h, and call_var_main.o depended on bam_utils.c rather than .h. The
# effect was that editing bam_utils.h left those objects stale, and because
# bam_chunk_t is defined there, a half-rebuilt tree segfaults on a struct
# layout mismatch instead of failing to compile.
-include $(OBJS:.o=.d)

.PHONY: all hts_all abpoa_all wfa2_all clean clean_all clean_hts clean_abpoa clean_wfa2

clean:
	rm -f $(SRC_DIR)/*.o $(SRC_DIR)/*.d $(EDLIB_DIR)/src/*.d $(BIN) .build-mode .build-version
clean_all:
	rm -f $(SRC_DIR)/*.o $(SRC_DIR)/*.d $(EDLIB_DIR)/src/*.d $(BIN) $(HTSLIB) $(ABPOA_LIB) $(WFA2_LIB)
clean_hts:
	rm -f $(HTSLIB)
clean_abpoa:
	rm -f $(ABPOA_LIB) $(ABPOA_DIR)/src/*.o
clean_wfa2:
	rm -f $(WFA2_LIB)

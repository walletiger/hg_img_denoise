
CC?=gcc
CFLAGS?=-Ofast -fPIC -fPIE
CFLAGS+=-Iinclude/

CXX?=g++
CXXFLAGS?=-Ofast -fPIC -fPIE
CXXFLAGS+=-Iinclude/


LOCAL_OBJ_FILES := \
	src/hg_guided_denoise2.o           \
	src/hg_integral_like_filter.o      \
	src/hg_local_variance_denoise_arm_neon.o      \
	src/hg_clahe_gray8.o       \
	src/hg_guided_denoise2_arm_neon.o       \
	src/hg_guided_denoise.o    \
	src/hg_local_variance_denoise.o       

.cc.o:
	$(CXX) -c $(CXXFLAGS) $*.cc -o $*.o

.c.o:
	$(CC) -c $(CFLAGS) $*.c -o $*.o

all: lib/libhg_denoise.a

lib/libhg_denoise.a: $(LOCAL_OBJ_FILES)
	$(AR) $(ARFLAGS) $@ $(LOCAL_OBJ_FILES)

test: bin/test_denoise

bin/test_denoise:
	$(CXX) $(CXXFLAGS) -Iinclude test/test_denoise.c lib/libhg_denoise.a  -o bin/test_denoise

clean:
	rm -fr src/*.o  rm -fr bin/* rm -fr lib/*




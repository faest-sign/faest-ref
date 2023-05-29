CC?=gcc
CXX?=g++
CFLAGS+=-g -O3 -std=c11
CXXFLAGS+=-g -O3 -std=c++17
CPPFLAGS+=-I. -Isha3 -DNDEBUG

SOURCES=$(filter-out randomness.c,$(wildcard *.c)) $(wildcard sha3/*.c) $(wildcard sha3/*.s)
EXECUTABLE_APITEST=tests/api_test
EXECUTABLE_NISTKAT=NIST-KATs/PQCgenKAT_sign
LIBFAEST=libfaest.a
NISTKAT_SOURCES=$(wildcard NIST-KATs/*.c) $(wildcard NIST-KATs/*.cpp)

all: $(LIBFAEST) $(EXECUTABLE_APITEST) $(EXECUTABLE_NISTKAT)

$(LIBFAEST): $(SOURCES:.c=.o) $(SOURCES:.s=.o)
	ar rcs $@ $^

%.c.o: %.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@

%.cpp.o: %.cpp
	$(CC) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

$(EXECUTABLE_NISTKAT): CPPFLAGS+=-DHAVE_RANDOMBYTES
$(EXECUTABLE_NISTKAT): CFLAGS+=-Wno-sign-compare -Wno-unused-but-set-variable -Wno-unused-parameter -Wno-unused-result

$(EXECUTABLE_APITEST): $(EXECUTABLE_APITEST).c.o $(LIBFAEST) randomness.c.o
	$(CC) $(CPPFLAGS) $(LDFLAGS) $^ -o $@

$(EXECUTABLE_NISTKAT): $(addsuffix .o, $(NISTKAT_SOURCES)) $(LIBFAEST) randomness.c.o
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $^ -lcrypto -o $@

clean:
	rm -f $(wildcard *.o) $(wildcard sha3/*.o) $(LIBFAEST) $(EXECUTABLE_APITEST) $(EXECUTABLE_NISTKAT)

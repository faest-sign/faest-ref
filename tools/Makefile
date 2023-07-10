CC?=gcc
CXX?=g++
CFLAGS+=-g -O3 -std=gnu11 -march=native
CPPFLAGS+=-I. -Isha3 -DNDEBUG -DHAVE_OPENSSL -DFAEST_EXPORT=

SOURCES=$(filter-out randomness.c,$(wildcard *.c)) $(wildcard sha3/*.c) $(wildcard sha3/*.s)
EXECUTABLE_APITEST=tests/api_test
EXECUTABLE_NISTKAT=NIST-KATs/PQCgenKAT_sign
LIBFAEST=libfaest.a
NISTKAT_SOURCES=$(wildcard NIST-KATs/*.c)

all: $(LIBFAEST) $(EXECUTABLE_APITEST) $(EXECUTABLE_NISTKAT)

$(LIBFAEST): $(SOURCES:.c=.o) $(SOURCES:.s=.o)
	ar rcs $@ $^

%.c.o: %.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@

$(EXECUTABLE_NISTKAT): CPPFLAGS+=-DHAVE_RANDOMBYTES
$(EXECUTABLE_NISTKAT): CFLAGS+=-Wno-sign-compare -Wno-unused-but-set-variable -Wno-unused-parameter -Wno-unused-result

$(EXECUTABLE_APITEST): $(EXECUTABLE_APITEST).c.o $(LIBFAEST) randomness.c.o
	$(CC) $(CPPFLAGS) $(LDFLAGS) $^ -lcrypto -o $@

$(EXECUTABLE_NISTKAT): $(addsuffix .o, $(NISTKAT_SOURCES)) $(LIBFAEST) randomness.c
	$(CC) $(CPPFLAGS) $(LDFLAGS) $^ -lcrypto -o $@

clean:
	rm -f \
		$(wildcard *.o) \
		$(wildcard sha3/*.o) \
		$(wildcard NIST-KATs/*.o) \
		$(wildcard tests/*.o) \
		$(LIBFAEST) $(EXECUTABLE_APITEST) $(EXECUTABLE_NISTKAT)

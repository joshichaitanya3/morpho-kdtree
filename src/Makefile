src = $(wildcard *.c) $(wildcard */*.c)
obj = $(src:.c=.dylib)

LDFLAGS  = -lm -lcblas -llapack -lcxsparse -L/opt/homebrew/lib
CFLAGS   =  -undefined dynamic_lookup -dynamiclib $(EXTCFLAGS) -O3 -I/opt/homebrew/include -I/usr/local/include/morpho

fastkdtree: 
	$(CC) -o $@.dylib $(src) $(LDFLAGS) $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(obj) 

.PHONY: install
install: fastkdtree
	cp fastkdtree.dylib ../lib/fastkdtree.dylib
	make clean 
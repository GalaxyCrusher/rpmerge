#Name of C compiler
CCMPL = gcc

rpmerge_v1.0: rpmerge_v1.0.o
	$(CCMPL) -o rpmerge_v1.0 rpmerge_v1.0.o -lm -g


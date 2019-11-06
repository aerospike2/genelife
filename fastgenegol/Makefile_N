
src = $(filter-out %sumtest.c %fastgenegol.c %fastgenelife.c %fastgol.c %fastgolstat.c %subgenelife_old.c %actgenelife.c , $(wildcard *.c))
obj = $(src:.c=.o)

LDFLAGS = -lGL

genelife_nopython: $(obj)
	$(CC) -o $@ $^ $(LDFLAGS)

foo:
	@echo ${src}


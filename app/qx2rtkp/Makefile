PROGNAME = qx2rtkp

LIBSRC = ../../src
OPTS    = -DTRACE -DENAGLO -DENAQZS -DENAGAL -DNFREQ=3

CC = gcc
LD = gcc

CFLAGS  = -Wall -ansi -pedantic -std=gnu99 -Wno-unused-but-set-variable -I$(LIBSRC) $(OPTS) -g -DTRACE
#CFLAGS += -I..
LDLIBS  = -lm -lrt


OBJDIR = build

OUT = $(OBJDIR)/$(PROGNAME)

RT_SRCS += rtkcmn.c rinex.c rtkpos.c postpos.c solution.c
RT_SRCS += lambda.c geoid.c sbas.c preceph.c pntpos.c ephemeris.c options.c
RT_SRCS += ppp.c ppp_ar.c rtcm.c rtcm2.c rtcm3.c rtcm3e.c ionex.c qzslex.c

SRCS= qx2rtkp.c $(addprefix $(LIBSRC)/,$(RT_SRCS))

OBJS=$(SRCS:.c=.o)


all: pre clean $(OUT)

$(OUT): pre $(OBJS)
	$(LD) $(LDLIBS) -o $(OUT) $(addprefix $(OBJDIR)/, $(notdir $(OBJS)))
	@echo "\n--- BUILD OK"
	
%.o : %.c
	$(CC) -c $(CFLAGS) $< -o $(OBJDIR)/$(notdir $@)

   
clean:
	rm -rf $(OBJDIR)/*
	@echo "CLEANED"
	
	@echo SRCS : $(SRCS)
	@echo OBJS : $(OBJS)
	@echo OBJS_DIRS : $(notdir $(OBJS))
	

pre:
	@mkdir -p $(OBJDIR)
    
    

.PHONY: pre clean all


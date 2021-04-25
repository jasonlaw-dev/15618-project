TARGET = main
OBJS += main.o
# OBJS += main.o

CC = mpic++
#CFLAGS = -Wall -Werror -DDEBUG -g # debug flags
CFLAGS = -lstdc++ -Wall -Werror -g -O3 # release flags
CFLAGS += -MMD -MP
LDFLAGS += $(LIBS)

default:	$(TARGET)
all:		$(TARGET)

$(TARGET):	$(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ ../FreeImage/Dist/libfreeimage.a $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) -I../FreeImage/Dist -c -o $@ $<

DEPS = $(OBJS:%.o=%.d)
-include $(DEPS)

clean:
	rm $(TARGET) $(OBJS) $(DEPS) || true

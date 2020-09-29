TARGET=helloworld
include ../shared/Makefile.shared

$(TARGET):$(TARGET).o
	$(CC) $(TARGET).o -lm $(SNIPER_LDFLAGS) -o $(TARGET)


run_$(TARGET):
	../../run-sniper -c gainestown --roi -- ./helloworld 


CLEAN_EXTRA=viz

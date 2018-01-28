TARGET = radar_basic

all: $(TARGET)

radar_basic:		main.o pulse_train.o range_resolution.o doppler_frequency.o
	$(CC) -o radar_basic main.o pulse_train.o range_resolution.o doppler_frequency.o -lm

pulse_train:		pulse_train.o
	$(CC) -c pulse_train.c

range_resolution:	range_resolution.o
	$(CC) -c range_resolution.c

doppler_frequency:	doppler_frequency.o
	$(CC) -c doppler_frequency.c

clean:
	rm -rf $(TARGET)
	rm -rf *.o

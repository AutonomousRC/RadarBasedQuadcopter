TARGET = radar_basic pulse_train range_resolution doppler_frequency

all: $(TARGET)

radar_basic:		main.o pulse_train.o range_resolution.o doppler_frequency.o
	$(CC) -o $@ main.o pulse_train.o range_resolution.o doppler_frequency.o

pulse_train:		pulse_train.o
	$(CC) -o $@ pulse_train.o

range_resolution:	range_resolution.o
	$(CC) -o $@ range_resolution.o

doppler_frequency:	doppler_frequency.o
	$(CC) -o $@ doppler_frequency.o

clean:
	rm -rf $(TARGET)
	rm -rf *.o

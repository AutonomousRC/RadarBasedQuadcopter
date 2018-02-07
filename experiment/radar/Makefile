TARGET = radar_basic snr_dr_rcs
SRCS = pulse_train.c range_resolution.c doppler_frequency.c radar_range_equation.c power_aperture.c self_screening_jammer.c sdjpn.c burn_through_range.c stand_off_jammer.c math_tech.c ogl_helper.c
OBJS = $(SRCS:%.c=%.o)

all: $(TARGET)

snr_dr_rcs:		snr_dr_rcs.o radar_range_equation.o math_tech.o ogl_helper.o
	$(CC) -o snr_dr_rcs snr_dr_rcs.o radar_range_equation.o math_tech.o ogl_helper.o -lm -lGL -lglut -lGLU

radar_basic:		main.o pulse_train.o range_resolution.o doppler_frequency.o radar_range_equation.o power_aperture.o self_screening_jammer.o sdjpn.o burn_through_range.o stand_off_jammer.o math_tech.o
	$(CC) -o radar_basic main.o math_tech.o pulse_train.o range_resolution.o doppler_frequency.o radar_range_equation.o power_aperture.o self_screening_jammer.o sdjpn.o burn_through_range.o stand_off_jammer.o -lm

### First make it simple

#pulse_train:		pulse_train.o
#	$(CC) -c pulse_train.c
#
#range_resolution:	range_resolution.o
#	$(CC) -c range_resolution.c
#
#doppler_frequency:	doppler_frequency.o
#	$(CC) -c doppler_frequency.c
#
#radar_range_equation:	radar_range_equation.o
#	$(CC) -c radar_range_equation.c
#
#power_aperture:		power_aperture.o
#	$(CC) -c power_aperture.c
#
#self_screening_jammer:	self_screening_jammer.o
#	$(CC) -c self_screening_jammer.c
#
#sdjpn:			sdjpn.o
#	$(CC) -c sdjpn.c
#
#burn_through_range:	burn_through_range.o
#	$(CC) -c burn_through_range.c
#
#stand_off_jammer:	stand_off_jammer.o
#	$(CC) -c stand_off_jammer.c
#
#math_tech:		math_tech.o
#	$(CC) -c math_tech.c
#
#ogl_helper:		ogl_helper.o
#	$(CC) -c ogl_helper.c

clean:
	rm -rf $(TARGET)
	rm -rf *.o

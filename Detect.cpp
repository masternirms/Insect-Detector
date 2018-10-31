#include "Detect.h"
#include <portaudio.h>
#include <sndfile.hh>
#include <fftw3.h>
#include <sndfile.h>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

static int callback_during_recording(const void *input_buffer, void *output_buffer, unsigned long frames_per_buffer,
	const PaStreamCallbackTimeInfo* time_info, PaStreamCallbackFlags status_flags, void *user_data)
{
	// Prevent unused variable warnings.
	(void)output_buffer;
	(void)time_info;
	(void)status_flags;

	recorded_data *data = (recorded_data *)user_data;
	unsigned long frames_left = data->max_frame_index - data->frame_index;
	int finished;
	long frames_to_calculate;
	if (frames_left < frames_per_buffer)
	{
       frames_to_calculate = frames_left;
       finished = paComplete;
	}
	else
	{
		frames_to_calculate = frames_per_buffer;
		finished = paContinue;
	}

	const SAMPLE *read_buffer = (const SAMPLE*)input_buffer;
	SAMPLE *write_buffer = &data->samples[data->frame_index * NUM_CHANNELS];
	if (input_buffer == NULL )
	{
		for (long i = 0; i < frames_to_calculate; i++)
		{
			*write_buffer++ = SAMPLE_SILENCE;  // Left
			if (NUM_CHANNELS == 2)
			{
				*write_buffer++ = SAMPLE_SILENCE;  // Right
			}
		}
	}
	else
	{
		for (long i = 0; i < frames_to_calculate; i++)
		{
			*write_buffer++ = *read_buffer++;  // Left
			if (NUM_CHANNELS == 2)
			{
				*write_buffer++ = *read_buffer++; // Right
			}
		}
	}

	data->frame_index += frames_to_calculate;

	return finished;
}

bool record_stream()
{
	recorded_data data;
	data.max_frame_index = RECORDING_DURATION * SAMPLE_RATE;
	data.frame_index = 0;
	int num_samples = data.max_frame_index * NUM_CHANNELS;
	int num_bytes = num_samples * sizeof(SAMPLE);
	data.samples = new SAMPLE[num_bytes](); // Initialize to zero.

	if( data.samples == NULL )
	{
		cout << "Could not allocate buffer for recording..." << endl;
		return false;
	}

	if (Pa_Initialize() != paNoError)
	{
		cout << "Failed to initialize port audio..." << endl;
		return false;
	}

	PaStreamParameters  inputParameters;
	inputParameters.device = Pa_GetDefaultInputDevice();
	if (inputParameters.device == paNoDevice)
	{
		cout << "Failed to open audio device..." << endl;
		return false;
	}

	inputParameters.channelCount = NUM_CHANNELS;
	inputParameters.sampleFormat = PA_SAMPLE_TYPE;
	inputParameters.suggestedLatency = Pa_GetDeviceInfo(inputParameters.device)->defaultLowInputLatency;
	inputParameters.hostApiSpecificStreamInfo = NULL;

	PaStream *stream;
	PaError error = Pa_OpenStream(&stream, &inputParameters, NULL, SAMPLE_RATE, FRAMES_PER_BUFFER, paClipOff,
		callback_during_recording, &data);
	if (error != paNoError)
	{
		cout << "Failed to open audio stream..." << endl;
		return false;
	};

	error = Pa_StartStream(stream);
	if (error != paNoError)
	{
		cout << "Failed to start audio stream..." << endl;
		return false;
	};

	cout << "Starting recording for " << RECORDING_DURATION << " seconds. Play back your mosquito sound..." << endl;

	while((error = Pa_IsStreamActive(stream)) == PA_STREAM_ACTIVE)
	{
		Pa_Sleep(ONE_SECOND_SLEEP);
	}
	if (error < paNoError)
	{
		cout << "Error in stream.." << endl;
		return false;
	}

	error = Pa_CloseStream(stream);
	if (error != paNoError)
	{
		cout << "Error in closing stream.." << endl;
		return false;
	}

	SndfileHandle sound_file = SndfileHandle ("recorded.wav", SFM_WRITE, SF_FORMAT_WAV | SF_FORMAT_PCM_16, NUM_CHANNELS,
		SAMPLE_RATE);
	sound_file.write (data.samples, num_samples) ;

	delete[] data.samples;
	Pa_Terminate();

	return true;
}

SndfileHandle open_file(const char *file_name)
{
	SndfileHandle sound_file = SndfileHandle(file_name);
	cout << "Failed to open sound file..." << endl;

	return sound_file;
}

sf_count_t read_file(SndfileHandle sound_file, double *buffer, int count)
{
	sf_count_t read_count = sound_file.read(buffer, count) ;
	if (read_count == count)
	{
		cout << "Read file..." << endl;
	}

	return read_count;
}

void do_fft(const double *input, fftw_complex *result, int count)
{
	fftw_plan plan = fftw_plan_dft_r2c_1d(count, (double *)input, result, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
}

void get_power_spectrum(double *amplitudes, double *frequencies, const fftw_complex *fft_result,
	int count, double bin_width, ofstream& csv_file)
{
	for (int i = 0; i < count; i++)
	{
		double normalized_bin_magazine = 2 * sqrt(fft_result[i][0] * fft_result[i][0] +
			fft_result[i][1] * fft_result[i][1]) / (count /* * 0.54*/); // 0.54 correction for a Hanning Window
		double amplitude = 20 * log10(normalized_bin_magazine);
		amplitudes[i] = amplitude;

		double frequency = i * bin_width;
		frequencies[i] = frequency;
		if (csv_file)
		{
			csv_file << i << "," << frequency << "," << amplitude << endl;
		}
	}
}

double get_fundamental_frequency_index(const double *amplitudes, int count, int harmonics)
{
	double *amplitude_sums = new double[count];
	amplitude_sums[0] = amplitudes[0];

	for (int i = 1; i <= count; i++)
	{
		double sum = 0;
		for (int j = 1; j <= harmonics; j++)
		{
			sum += amplitudes[i * j];
		}

		amplitude_sums[i] = sum;
	}

	int max = 1;
	double max_sum = amplitude_sums[1];
	for (int i = 1; i <= count; i++)
	{
		if (amplitude_sums[i] > max_sum)
		{
			max = i;
			max_sum = amplitude_sums[i];
		}
	}

	delete [] amplitude_sums;

	return max;
}

void get_amplitudes_of_interest(int index, const double *amplitudes, double *harmonic_amplitudes, int harmonics)
{
	for (int i = 1; i <= harmonics; i++)
	{
		harmonic_amplitudes[i - 1] = amplitudes[index * i];
	}

	for (int i = 0; i <= harmonics; i++)
	{
		int start = index / 2;
		harmonic_amplitudes[harmonics + i] = amplitudes[start + index * i];
	}
}

bool determine_mosquito_presence(const double *harmonic_amplitudes, int harmonics)
{
	// Difference between 1st harmonic and adjacent inter-peaks: lower threshold - 26dB.
	double harmonic_value = harmonic_amplitudes[0];
	double lower_inter_peak_value = harmonic_amplitudes[harmonics];
	//cout << harmonic_value - lower_inter_peak_value << endl;
	if ((harmonic_value - lower_inter_peak_value) < 20 /*26*/)
	{
		return false;
	}
	double upper_inter_peak_value = harmonic_amplitudes[harmonics + 1];
	//cout << harmonic_value - upper_inter_peak_value << endl;
	if ((harmonic_value - upper_inter_peak_value) < 20 /*26*/)
	{
		return false;
	}

	// Difference between 2nd harmonic and adjacent inter-peaks: lower threshold - 8dB.
	harmonic_value = harmonic_amplitudes[1];
	lower_inter_peak_value = harmonic_amplitudes[harmonics + 1];
	// cout << harmonic_value - lower_inter_peak_value << endl;
	if ((harmonic_value - lower_inter_peak_value) < 8)
	{
		return false;
	}
	upper_inter_peak_value = harmonic_amplitudes[harmonics + 2];
	// cout << harmonic_value - upper_inter_peak_value << endl;
	if ((harmonic_value - upper_inter_peak_value) < 8)
	{
		return false;
	}

	// Difference between 3rd harmonic and adjacent inter-peaks: lower threshold - 4dB.
	harmonic_value = harmonic_amplitudes[2];
	lower_inter_peak_value = harmonic_amplitudes[harmonics + 2];
	//cout << harmonic_value - lower_inter_peak_value << endl;
	if ((harmonic_value - lower_inter_peak_value) < 4)
	{
		return false;
	}
	upper_inter_peak_value = harmonic_amplitudes[harmonics + 3];
	//cout << harmonic_value - upper_inter_peak_value << endl;
	if ((harmonic_value - upper_inter_peak_value) < 4)
	{
		return false;
	}

	// Difference between 4th harmonic and adjacent inter-peaks: lower threshold - 2dB.
	harmonic_value = harmonic_amplitudes[3];
	lower_inter_peak_value = harmonic_amplitudes[harmonics + 3];
	// cout << harmonic_value - lower_inter_peak_value << endl;
	if ((harmonic_value - lower_inter_peak_value) < 2)
	{
		return false;
	}
	upper_inter_peak_value = harmonic_amplitudes[harmonics + 4];
	// cout << harmonic_value - upper_inter_peak_value << endl;
	if ((harmonic_value - upper_inter_peak_value) < 2)
	{
		return false;
	}

	// Difference between 1st harmonic and 2nd harmonic: lower threshold - (-6)dB, upper threshold - 16dB.
	double harmonic_value1 = harmonic_amplitudes[0];
	double harmonic_value2 = harmonic_amplitudes[1];
	double difference = harmonic_value1 - harmonic_value2;
	if (difference < -6)
	{
		return false;
	}
	if (difference > 16)
	{
		return false;
	}

	// Difference between 2nd harmonic and 3rd harmonic: lower threshold - 4dB, upper threshold - 26dB.
	harmonic_value1 = harmonic_amplitudes[1];
	harmonic_value2 = harmonic_amplitudes[2];
	difference = harmonic_value1 - harmonic_value2;
	if (difference < 4)
	{
		return false;
	}
	if (difference > 26)
	{
		return false;
	}

	// Difference between 3rd harmonic and 4th harmonic: lower threshold - 8dB, upper threshold - 34dB.
	harmonic_value1 = harmonic_amplitudes[2];
	harmonic_value2 = harmonic_amplitudes[3];
	difference = harmonic_value1 - harmonic_value2;
	if (difference < 8)
	{
#if 0
		return false;
#endif
	}
	if (difference > 34)
	{
		return false;
	}

	return true;
}


int main(void)
{
	if (!record_stream())
	{
		return 1;
	}
	SndfileHandle sound_file =  open_file("recorded.wav");; // RAII takes care of destroying this object.
	if (!sound_file)
	{
		return 1;
	}

	ofstream csv_file;
#ifdef DEBUG
	csv_file.open("sound_data.csv", ios::out);
	if (!csv_file)
	{
		cout << "Failed to create CSV file..." << endl;
		return 1;
	}
#endif // DEBUG

	while (true)
	{
		double samples[SAMPLE_COUNT];
		if (read_file(sound_file, samples, SAMPLE_COUNT) < SAMPLE_COUNT)
		{
			cout << "Processing complete..." << endl;
			return 0;
		}

		fftw_complex fft_result[SAMPLE_COUNT];
		do_fft(samples, fft_result, SAMPLE_COUNT);

		double amplitudes[SAMPLE_COUNT];
		double frequencies[SAMPLE_COUNT];
		double bin_width = sound_file.samplerate()/ SAMPLE_COUNT;
		get_power_spectrum(amplitudes, frequencies, fft_result, SAMPLES_OF_INTEREST, bin_width, csv_file);

		int fundamental_frequency_index = get_fundamental_frequency_index(amplitudes, BIN_COUNT_OF_INTEREST, HARMONICS_OF_INTEREST);
		cout << "Fundamental Frequency: " << frequencies[fundamental_frequency_index] << endl;
		if (fundamental_frequency_index > 1)
		{
			double harmonic_amplitudes[HARMONICS_OF_INTEREST * 2 + 1];
			get_amplitudes_of_interest(fundamental_frequency_index, amplitudes, harmonic_amplitudes, HARMONICS_OF_INTEREST);

			bool result = determine_mosquito_presence(harmonic_amplitudes, HARMONICS_OF_INTEREST);
			cout << "Mosquito " << (result ? "is " : "not ") << "detected..." << endl;
		}
	}
}

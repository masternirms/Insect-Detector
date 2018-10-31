#define NUM_CHANNELS 1 // Mono
#define SAMPLE_RATE 48000
#define FRAMES_PER_BUFFER 512
#define PA_SAMPLE_TYPE  paFloat32
#define RECORDING_DURATION 30 // 30 seconds 
#define PA_STREAM_ACTIVE 1
#define ONE_SECOND_SLEEP 1000
#define SAMPLE_SILENCE 0.0f
#define SAMPLE_COUNT 4096
#define SAMPLES_OF_INTEREST 683 /* AC97 sampling is at 48kHz; not 44.1kHz; 683 samples correspond to 8kHz */
#define BIN_COUNT_OF_INTEREST 148
#define HARMONICS_OF_INTEREST 4

typedef float SAMPLE;

typedef struct
{
   unsigned long frame_index;  // Index into the frames array
   unsigned long max_frame_index;
   SAMPLE *samples;
} recorded_data;

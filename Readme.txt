sudo apt-get install libfftw3-dev
sudo apt-get install libsndfile-dev
sudo apt-get install portaudio19-dev

gcc detect.cpp -lstdc++ -lm -lfftw3 -lsndfile -lportaudio -Wall -o detect
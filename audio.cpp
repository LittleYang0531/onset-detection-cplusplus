#include<bits/stdc++.h>
#include"libaudio.h"
using namespace std;
int main(int argc, char** argv) {
    wavFile wav;
    wav.open(argv[1]);
    wav.outputBasicInfo();
    audioData audio = audioData(wav, true, 20);
    specturmData sp = specturmData(audio, 48, true, 20);
    beatData beat = beatData(sp, 10, 1.1);
    beat.drawGraph("test.svg");
    beat.onsetDetection(0.1);
    int beatNum = beat.fetchBeatNum();
    float* beatTime = beat.fetchBeatTime();
    cout << "Beat Number: " << beatNum << endl;

    ofstream fout(argv[2]);
	stringstream notes;
	for (int i = 0; i < beatNum; i++) {
		notes << "{\"type\": \"Single\", \"lane\": " << rand() % 7 
              << ", \"beat\": " << 2 * beatTime[i] << "}" << (i != beatNum - 1 ? ", " : "");
	}
	fout << "[{\"type\": \"System\", \"cmd\": \"BPM\", \"beat\": 0, \"bpm\": 90}, ";
	fout << notes.str() << "]";
}
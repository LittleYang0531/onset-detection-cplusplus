#include<bits/stdc++.h>
using namespace std;
const float PI = acos(-1);

struct comp {
	float real = 0, virt = 0;
	comp(){};
	comp(float real, float virt): real(real), virt(virt){};
	comp operator + (const comp& a) const {
		return comp(real + a.real, virt + a.virt);
	}
	comp operator - (const comp& a) const {
		return comp(real - a.real, virt - a.virt);
	}
	comp operator * (const comp& a) const {
		return comp(real * a.real - virt * a.virt, real * a.virt + virt * a.real);
	}
};

comp* complexize(float* st, float* en) {
	comp* res = new comp[en - st];
	float* orig = st; comp* goal = res;
	for (; orig != en; orig++, goal++) *goal = comp(*orig, 0);
	return res;
}

float* floatize(comp* st, comp* en) {
	float* res = new float[en - st];
	comp* orig = st; float* goal = res;
	for (; orig != en; orig++, goal++) *goal = sqrt(orig->real * orig->real + orig->virt * orig->virt);
	return res;
}

comp getw(int n, int k) {
	return comp(cos(2 * PI * k / n), sin(2 * PI * k / n));
}

void fft(comp* st, comp* en, bool inv = false) {
	if (st == en) return;
	int len = en - st;
	assert((len & (-len)) == len);
	comp* odd = new comp[len / 2];
	comp* even = new comp[len / 2];
	for (int i = 0; i < len; i++) {
		if (i % 2 == 0) odd[i / 2] = *(st + i);
		else even[i / 2] = *(st + i);
	}
	fft(odd, odd + len / 2, inv);
	fft(even, even + len / 2, inv);
	comp w = comp(1, 0);
	for (int i = 0; i < len / 2; i++) {
		*(st + i) = odd[i] + w * even[i];
		*(st + i + len / 2) = odd[i] - w * even[i];
		w = w * (inv ? getw(-1, len) : getw(1, len));
	}
	free(odd), free(even);
}

struct wavfile {
	string ChunkID; int ChunkSize; string Format; 
	string SubChunk1ID; int SubChunk1Size;
	int AudioFormat; int NumChannels; int SampleRate;
	int ByteRate; int BlockAlign; int BitsPerSample;
	struct SubChunk {
		string ID;
		int Size;
		char* Data;
	} SubChunk[1024];
	int SubChunkSize = 1;
	
	void open(const char* path) {
		// 文件读取
		ifstream fin(path);
		fin.seekg(0, ios::end);
		int length = fin.tellg();
		fin.seekg(0, ios::beg);
		char* tmp = new char[length];
		fin.read(tmp, length);
		int base;

		// 写入 ChunkID
		for (int i = 0; i < 4; i++) ChunkID += tmp[i];
		// 写入 ChunkSize
		base = 1;
		for (int i = 4; i < 8; i++) ChunkSize += base * (unsigned char)tmp[i], base <<= 8;
		// 写入 Format
		for (int i = 8; i < 12; i++) Format += tmp[i];
		// 写入 SubChunk1ID
		for (int i = 12; i < 16; i++) SubChunk1ID += tmp[i];
		// 写入 SubChunk1Size
		base = 1;
		for (int i = 16; i < 20; i++) SubChunk1Size += base * (unsigned char)tmp[i], base <<= 8;
		// 写入 AudioFormat
		base = 1;
		for (int i = 20; i < 22; i++) AudioFormat += base * (unsigned char)tmp[i], base <<= 8;
		// 写入 NumChannels
		base = 1;
		for (int i = 22; i < 24; i++) NumChannels += base * (unsigned char)tmp[i], base <<= 8;
		cout << "Channel Number: " << NumChannels << endl;
		// 写入 SampleRate
		base = 1;
		for (int i = 24; i < 28; i++) SampleRate += base * (unsigned char)tmp[i], base <<= 8;
		cout << "Sample Rate: " << SampleRate << endl;
		// 写入 ByteRate
		base = 1;
		for (int i = 28; i < 32; i++) ByteRate += base * (unsigned char)tmp[i], base <<= 8;
		cout << "Byte Rate: " << ByteRate << endl;
		// 写入 BlockAlign
		base = 1;
		for (int i = 32; i < 34; i++) BlockAlign += base * (unsigned char)tmp[i], base <<= 8;
		// 写入 BitsPerSample
		base = 1;
		for (int i = 34; i < 36; i++) BitsPerSample += base * (unsigned char)tmp[i], base <<= 8;
		cout << "Bits per Sample: " << BitsPerSample << endl;
		int summary = ChunkSize, pt = 36; summary -= 36;
		while (summary > 0) {
			SubChunkSize++;
			// 写入 SubChunkID
			for (int i = pt; i < pt + 4; i++) SubChunk[SubChunkSize].ID += tmp[i];
			// 写入 SubChunkSize
			pt += 4, base = 1;
			for (int i = pt; i < pt + 4; i++) SubChunk[SubChunkSize].Size += base * (unsigned char)tmp[i], base <<= 8;
			cout << "Sub Chunk " << SubChunk[SubChunkSize].ID << " Size: " << SubChunk[SubChunkSize].Size << endl;
			// 写入 Data 
			pt += 4; 
			SubChunk[SubChunkSize].Data = new char[SubChunk[SubChunkSize].Size];
			memcpy(SubChunk[SubChunkSize].Data, tmp + pt, SubChunk[SubChunkSize].Size);
			pt += SubChunk[SubChunkSize].Size, summary -= SubChunk[SubChunkSize].Size + 8;
		}
	}

	float* leftChannel;
	float* rightChannel;
	int SampleSize;
	void waveChannel() {
		int perSampleSize = NumChannels * BitsPerSample / 8;
		SampleSize = SubChunk[3].Size / perSampleSize;
		leftChannel = new float[SampleSize];
		rightChannel = new float[SampleSize];
		int pt = 0;
		for (int i = 0; i < SubChunk[3].Size; i += perSampleSize) {
			int base = 1, leftSample = 0, rightSample = 0;
			// 计算左右声道样本数据
			for (int j = 0; j < BitsPerSample / 8; j++)
				leftSample += base * (unsigned char)SubChunk[3].Data[i + j], base <<= 8;
			base = 1;
			for (int j = 0; j < BitsPerSample / 8; j++)
				rightSample += base * (unsigned char)SubChunk[3].Data[i + BitsPerSample / 8 + j], base <<= 8;
			
			// 对样本数据取补
			if (leftSample & (1 << (BitsPerSample - 1))) leftSample = -((~(leftSample - 1)) & ((1 << BitsPerSample) - 1));
			if (rightSample & (1 << (BitsPerSample - 1))) rightSample = -((~(rightSample - 1)) & ((1 << BitsPerSample) - 1));

			// 将样本数据转化为浮点数据
			leftChannel[pt] = 1.0 * leftSample / ((1 << BitsPerSample) - 1);
			rightChannel[pt++] = 1.0 * rightSample / ((1 << BitsPerSample) - 1);
		}
	}

	int windowSize;
	int windowNum;

	// 可调参数
	int averageSize = 20;
	float multiplier = 1;
	float limit = 0.5;
	// 可调参数

	float* prunnedSpecturm;
	void getSpectrum(int goalRate) {
		cout << "Goal Rate: " << goalRate << endl;
		windowSize = round(1.0 * SampleRate / goalRate);
		cout << "Sample Number per Window: " << windowSize << endl;
		windowNum = ceil(1.0 * SampleSize / windowSize);
		cout << "Window Number: " << windowNum << endl;
		cout << "Average Size: " << averageSize << endl;
		cout << "Multiply Times: " << multiplier << endl;
		int fftSize = windowSize, processNum = 0, percentNow = 0;
		while ((fftSize & (-fftSize)) != fftSize) fftSize++;
		float* lastSpectrum = NULL; float* Specturm; float* SF;
		SF = new float[windowNum]; float specturmMin = 1e18;
		for (int i = 0; i < SampleSize; i += windowSize) {
			float* tmp = new float[fftSize];
			for (int j = 0; j < fftSize; j++) tmp[j] = 0;
			for (int j = 0; j < windowSize && i + j < SampleSize; j++)
				tmp[j] = (leftChannel[i + j] + rightChannel[i + j]) / NumChannels; // 声道合并
			comp* complex = complexize(tmp, tmp + fftSize);
			fft(complex, complex + fftSize);
			Specturm = floatize(complex, complex + fftSize);
			float newSF = 0;
			for (int j = 0; j < fftSize; j++) {
				float value = Specturm[j] - (i != 0 ? lastSpectrum[j] : 0);
				if (value >= 0) newSF += value; 
			}
			free(tmp); free(complex); free(lastSpectrum);
			lastSpectrum = Specturm;
			SF[processNum] = newSF; processNum++;
			specturmMin = min(specturmMin, SF[processNum]);
			if (processNum * 10 / windowNum != percentNow) {
				percentNow = processNum * 10 / windowNum;
				cout << "Processed: " << percentNow * 10 << "%" << endl;
			}
		}
		float* threshold = new float[windowNum];
		float thresholdMin = 1e18;
		for (int i = 0; i < windowNum; i++) {
			int leftBorder = max(0, i - averageSize / 2);
			int rightBorder = min(windowNum - 1, i + averageSize / 2);
			float summarySF = 0;
			for (int j = leftBorder; j < rightBorder; j++) summarySF += SF[j];
			summarySF /= rightBorder - leftBorder;
			threshold[i] = summarySF * multiplier;
			thresholdMin = min(thresholdMin, threshold[i]);
		}
		prunnedSpecturm = new float[windowNum];
		for (int i = 0; i < windowNum; i++) {
			if (threshold[i] <= SF[i]) prunnedSpecturm[i] = SF[i] - threshold[i];
			else prunnedSpecturm[i] = 0;
		}
	}

	int* beats;
	int beatSize = 0;
	void beatDetection() {
		float maxSpecturm = 0, whereSpecturm;
		for (int i = 0; i < windowNum; i++) {
			if (prunnedSpecturm[i]) {
				if (prunnedSpecturm[i] > maxSpecturm) whereSpecturm = i;
				maxSpecturm = max(maxSpecturm, prunnedSpecturm[i]);
			} else {
				if (maxSpecturm > limit) beatSize++;
				maxSpecturm = 0;
			}
		}
		int solvedNum = 0;
		cout << "Beats Number: " << beatSize << endl;	
		beats = new int[beatSize];
		maxSpecturm = 0;
		for (int i = 0; i < windowNum; i++) {
			if (prunnedSpecturm[i]) {
				if (prunnedSpecturm[i] > maxSpecturm) whereSpecturm = i;
				maxSpecturm = max(maxSpecturm, prunnedSpecturm[i]);
			} else {
				if (maxSpecturm > limit) beats[solvedNum++] = whereSpecturm;
				maxSpecturm = 0;
			}
		}
	}
};

int main(int argc, char** argv) {
	if (argc < 3) {
		cout << "Usage: " << argv[0] << " [wav] [json]" << endl;
		return 3;
	}
	wavfile wav;
	wav.open(argv[1]);
	wav.waveChannel();
	wav.getSpectrum(72);
	wav.beatDetection();
	ofstream fout(argv[2]);
	
	int bpb = 4, timepoint = time(NULL);
	
	stringstream notes;
	for (int i = 0; i < wav.beatSize; i++) {
		notes << "{\"type\": \"Single\", \"lane\": " << rand() % 7 << ", \"beat\": " << 1.0 * wav.beats[i] / 36 << "}" << (i != wav.beatSize - 1 ? ", " : "");
	}
	fout << "[{\"type\": \"System\", \"cmd\": \"BPM\", \"beat\": 0, \"bpm\": 90}, ";
	fout << notes.str() << "]";
}

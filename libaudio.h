#ifndef LIBAUDIO_H
#define LIBAUDIO_H
#define int long long
#include<bits/stdc++.h>
using namespace std;

class pointer {
    private:
    int _pt = 0, _ptLimit = 0;
    string _ptMsg = "";

    public:
    int x = 0;
    pointer(int beg, int end, string msg = "") {
        init(beg, end, msg);
    }
    pointer(){}
    void init(int beg, int end, string msg = "") {
        _pt = beg, _ptLimit = end, _ptMsg = msg;
        if (beg > end) {
            cerr << "Failed to initialize pointer: "
                 << "begin address is greater than end address!" << endl;
            abort();
        }
        x = beg;
    }
    void move(int step = 1) {
        _pt += step;
        if (_pt > _ptLimit) {
            cerr << _ptMsg << endl;
            abort();
        }
        x = _pt;
    }
    pointer operator ++ () {
        move(1);
        return *this;
    }
    #undef int
    pointer operator ++ (int) {
        pointer origin = *this;
        move(1);
        return origin;
    }
    #define int long long
    void operator += (int x) {
        move(x);
    }
    operator int() {
        return x;
    }
}; 

class wavFile {
    private:
    bool _open = false;
    string _path = "";

    public:
    string chunkId; int chunkSize; string format;
    int audioFormat; int numChannels; int sampleRate;
    int byteRate; int blockAlign; int bitsPerSample;
    int subChunkSize = 0;
    struct subChunk {
        string id;
        int size;
        char* data;
    } subChunk[1024];

    wavFile(){}
    wavFile(const char* dest) {
        open(dest);
    }
    void open(const char* dest) {
        if (_open) {
            cerr << "Failed to open wav file: Have opened an audio file" << endl;
            abort();
        }
        
        // 文件读入
        ifstream fin(dest);
        if (!fin) {
            cerr << "Failed to open wav file: no such file or directory" << endl;
            abort();
        }
		fin.seekg(0, ios::end);
		int length = fin.tellg();
		fin.seekg(0, ios::beg);
		char* tmp = new char[length];
		fin.read(tmp, length);

        int base = 1; pointer pt(0, length, "Invalid wav file: chunk length is not enough");
        // 写入 chunkId
		for (int i = 0; i < 4; i++) chunkId += tmp[pt], pt++;
        if (chunkId != "RIFF") {
            cerr << "Invalid wav file: unknown chunkId" << endl;
            abort();
        }
		// 写入 chunkSize
		base = 1;
		for (int i = 0; i < 4; i++) chunkSize += base * (unsigned char)tmp[pt], base <<= 8, pt++;
        // 写入 format
		for (int i = 0; i < 4; i++) format += tmp[pt], pt++;
		int summary = chunkSize; summary -= 4;
        while (summary > 0) {
            subChunkSize++;
            // 写入 subChunkId
            for (int i = 0; i < 4; i++) subChunk[subChunkSize].id += tmp[pt], pt++;
            // 写入 subChunkSize
            base = 1; for (int i = 0; i < 4; i++) subChunk[subChunkSize].size += base * (unsigned char)tmp[pt], base <<= 8, pt++;
            // 写入 subChunkData
            subChunk[subChunkSize].data = new char[subChunk[subChunkSize].size];
            memcpy(subChunk[subChunkSize].data, tmp + pt, subChunk[subChunkSize].size);
            pt += subChunk[subChunkSize].size, summary -= subChunk[subChunkSize].size + 8;
        }
        _open = true;
        _path = dest;

        // 获取音频基本属性
        freshBasicInfo();
    }
    bool isOpen() {
        return _open;
    }

    void freshBasicInfo() {
        if (!_open) {
            cerr << "Failed to get basic info: no wav file was opened" << endl;
            abort();
        }

        // 查找 fmt 子 chunk
        int id = -1;
        for (int i = 1; i <= subChunkSize; i++) if (subChunk[i].id == "fmt ") id = i;
        if (id == -1) {
            cerr << "Failed to get basic info: no subchunk named fmt was found" << endl;
            abort();
        }

        int base = 1; pointer pt(0, subChunk[id].size, "Invalid wav file: subchunk length is not enough");
        // 写入 AudioFormat
		base = 1;
		for (int i = 0; i < 2; i++) audioFormat += base * (unsigned char)subChunk[id].data[pt], base <<= 8, pt++;
		// 写入 NumChannels
		base = 1;
		for (int i = 0; i < 2; i++) numChannels += base * (unsigned char)subChunk[id].data[pt], base <<= 8, pt++;
		// 写入 SampleRate
		base = 1;
		for (int i = 0; i < 4; i++) sampleRate += base * (unsigned char)subChunk[id].data[pt], base <<= 8, pt++;
		// 写入 ByteRate
		base = 1;
		for (int i = 0; i < 4; i++) byteRate += base * (unsigned char)subChunk[id].data[pt], base <<= 8, pt++;
		// 写入 BlockAlign
		base = 1;
		for (int i = 0; i < 2; i++) blockAlign += base * (unsigned char)subChunk[id].data[pt], base <<= 8, pt++;
		// 写入 BitsPerSample
		base = 1;
		for (int i = 0; i < 2; i++) bitsPerSample += base * (unsigned char)subChunk[id].data[pt], base <<= 8, pt++;
    }
    void outputBasicInfo(ostream& out = cout) {
        if (!_open) {
            cerr << "Failed to output basic info: no wav file was opened" << endl;
            abort();
        }

        out << "Audio info for \"" << _path << "\": " << endl;
        out << "================================" << endl;
        out << "Main Chunk Size = " << chunkSize << endl;
        out << "File Format = " << format << endl;
        out << "Audio Format = " << audioFormat << endl;
        out << "Channel Number = " << numChannels << endl;
        out << "Sample Rate = " << sampleRate << endl;
        out << "Byte Rate = " << byteRate << endl;
        out << "Sample Size = " << blockAlign << endl;
        out << "Channel Sample = " << bitsPerSample << endl;
        out << "================================" << endl;
    }
    string getPath() {
        return _path;
    }

    void close() {
        if (!_open) {
            cerr << "Failed to close wav file: no wav file was opened" << endl;
            abort();
        }

        _path = "";
        for (int i = 1; i <= subChunkSize; i++) {
            free(subChunk[i].data);
            subChunk[i].data = NULL;
        }
        chunkId = ""; chunkSize = 0; format = "";
        audioFormat = 0; numChannels = 0; sampleRate = 0;
        byteRate = 0; blockAlign = 0; bitsPerSample = 0;
        subChunkSize = 0;
        _open = false;
    }
};

class mp3File {

};

#define LEFT_CHANNEL false
#define RIGHT_CHANNEL true
class audioData {
    private:
    string _path = "";
    float* leftChannel;
    float* rightChannel;
    bool _full = false;

    public:
    int numChannels = 0;
    int sampleNum = 0;
    int sampleRate = 0;
    audioData(){}
    audioData(int numChannels, int sampleNum, float* leftChannel, float* rightChannel, int sampleRate, string _path = ""): 
        _full(true), numChannels(numChannels), sampleNum(sampleNum), 
        _path(_path), leftChannel(leftChannel), rightChannel(rightChannel),
        sampleRate(sampleRate) {};
    audioData(int numChannels, int sampleNum, int sampleRate, string _path = ""): 
        numChannels(numChannels), sampleNum(sampleNum), sampleRate(sampleRate), _path(_path), _full(true){
        leftChannel = new float[sampleNum];
        rightChannel = new float[sampleNum];
        memset(leftChannel, 0, sampleNum);
        memset(rightChannel, 0, sampleNum);
    }
    audioData(wavFile wav, bool processBar = true, int partNum = 20) {
        if (!wav.isOpen()) {
            cerr << "Failed to convert to audio data: no wav file was opened" << endl;
            abort();
        }

        // 查找 data 子 chunk
        int id = -1;
		for (int i = 1; i <= wav.subChunkSize; i++) if (wav.subChunk[i].id == "data") id = i;
        if (id == -1) {
            cerr << "Failed to convert to audio data: no subchunk named data was found" << endl;
            abort();
        }

        // 处理需要的数据
		sampleNum = wav.subChunk[id].size / wav.blockAlign;
		leftChannel = new float[sampleNum];
		rightChannel = new float[sampleNum];
        numChannels = wav.numChannels;
        sampleRate = wav.sampleRate;
        _path = wav.getPath();
        int percentNow = 0;

        // 绘制加载进度条
        if (processBar) {
            cout << "Converting to Audio Wave Data: ";
            cout << "["; for (int i = 0; i < partNum; i++) cout << "-"; cout << "] 0/" << partNum;
        }
        cout.flush();

		pointer pt(0, wav.subChunk[id].size, "\nFailed to convert to audio data: subchunk length is not enough");
        int samplePt = 0;
		for (int i = 0; i < wav.subChunk[id].size; i += wav.blockAlign) {
			int base = 1, leftSample = 0, rightSample = 0;
			// 计算左右声道样本数据
			for (int j = 0; j < wav.bitsPerSample / 8; j++)
				leftSample += base * (unsigned char)wav.subChunk[id].data[pt], base <<= 8, pt++;
			base = 1;
			for (int j = 0; j < wav.bitsPerSample / 8; j++)
				rightSample += base * (unsigned char)wav.subChunk[id].data[pt], base <<= 8, pt++;
			
			// 对样本数据取补
			if (leftSample & (1ll << (wav.bitsPerSample - 1))) leftSample = -((~(leftSample - 1)) & ((1ll << wav.bitsPerSample) - 1));
			if (rightSample & (1ll << (wav.bitsPerSample - 1))) rightSample = -((~(rightSample - 1)) & ((1ll << wav.bitsPerSample) - 1));

			// 将样本数据转化为浮点数据
			leftChannel[samplePt] = 1.0 * leftSample / double(1ll << (wav.bitsPerSample - 1));
			rightChannel[samplePt++] = 1.0 * rightSample / double(1ll << (wav.bitsPerSample - 1));

            if (processBar && samplePt * partNum / sampleNum != percentNow) {
                // 绘制进度条
                int siz = partNum - percentNow + 3 + to_string(percentNow).size() + to_string(partNum).size();
                for (int j = 0; j < siz; j++) cout << "\b";
				percentNow = samplePt * partNum / sampleNum;

                cout << "#";
                for (int j = 0; j < partNum - percentNow; j++) cout << "-";
                cout << "] " << percentNow << "/" << partNum;
                cout.flush();
			}
		}

        _full = true;
        if (processBar) cout << endl;
        cout.flush();
    }

    int size() {
        return sampleNum;
    }
    void modifyChannel(int sampleId, float leftChannelData, float rightChannelData) {
        if (sampleId < 0 || sampleId >= sampleNum) {
            cerr << "Failed to modify channel data: out of range" << endl;
            abort();
        }
        *(leftChannel + sampleId) = leftChannelData;
        *(rightChannel + sampleId) = rightChannelData;
    }
    void modifyChannel(bool channelId, int sampleId, float channelData) {
        if (sampleId < 0 || sampleId >= sampleNum) {
            cerr << "Failed to modify channel data: out of range" << endl;
            abort();
        }
        if (channelId == LEFT_CHANNEL) *(leftChannel + sampleId) = channelData;
        if (channelId == RIGHT_CHANNEL) *(rightChannel + sampleId) = channelData;
    }
    void updateSize(int goalSize) {
        float* newLeftChannel = new float[goalSize];
        float* newRightChannel = new float[goalSize];
        memset(newLeftChannel, 0, goalSize);
        memset(newRightChannel, 0, goalSize);
        memcpy(newLeftChannel, leftChannel, min(goalSize, sampleNum));
        memcpy(newRightChannel, rightChannel, min(goalSize, sampleNum));
        free(leftChannel), free(rightChannel);
        leftChannel = newLeftChannel, rightChannel = newRightChannel;
    }

    pair<float*, float*> fetchSampleData() {
        return {leftChannel, rightChannel};
    }
    pair<float, float> fetchSampleData(int sampleId) {
        if (sampleId < 0 || sampleId >= sampleNum) {
            cerr << "Failed to fetch channel data: out of range" << endl;
            abort();
        }
        return {*(leftChannel + sampleId), *(rightChannel + sampleId)};
    }

    bool isFull() {
        return _full;
    }

    void drawGraph(const char* dest, int goalRate = 72, float multiRate = 100) {
		int windowSize = round(1.0 * sampleRate / goalRate);
		int windowNum = ceil(1.0 * sampleNum / windowSize);
        ofstream fout(dest);
        if (!fout) {
            cerr << "Failed to draw wave graph: couldn't open file" << endl;
            abort();
        }

        // 计算窗口值
        float minChannel = 1e18;
        float* minValue = new float[windowNum];
        float* maxValue = new float[windowNum];
        for (int i = 0; i < sampleNum; i += windowSize) {
            float minSample = 1e18, maxSample = -1e18;
            for (int j = 0; j < windowSize && i + j < sampleNum; j++) {
                minSample = min(minSample, min(leftChannel[i + j], rightChannel[i + j]));
                maxSample = max(maxSample, max(leftChannel[i + j], rightChannel[i + j]));
            }
            minValue[i / windowSize] = minSample;
            maxValue[i / windowSize] = maxSample;
            minChannel = min(minChannel, minSample);
        }

        // 输出图像
		fout << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">" << endl;
        for (int i = 1; i < windowNum; i++) {
            if (i % 2 == 1) fout << "<line x1=\"" << i - 1 << "\" y1=\"" << (minValue[i - 1] - minChannel) * multiRate
                                 << "\" x2=\"" << i << "\" y2=\"" << (maxValue[i] - minChannel) * multiRate
                                 << "\" style=\"stroke: rgb(255, 0, 0);\"/>" << endl;
            else fout << "<line x1=\"" << i - 1 << "\" y1=\"" << (maxValue[i - 1] - minChannel) * multiRate
                      << "\" x2=\"" << i << "\" y2=\"" << (minValue[i] - minChannel) * multiRate
                      << "\" style=\"stroke: rgb(255, 0, 0);\"/>" << endl;
        }
		fout << "</svg>" << endl;
        fout.close();
    }
};

class specturmData {
    private:
    float PI = acos(-1);
    bool _full = false;
    float* spData;

    // fft 实现部分
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

    public:
    int rate = 0;
    int windowSize = 0;
    int windowNum = 0;

    specturmData(){};
    specturmData(audioData audio, int goalRate = 72, bool processBar = true, int partNum = 20) {
        if (!audio.isFull()) {
            cerr << "Failed to convert to specturm data: no audio data was filled" << endl;
            abort();
        }

        // 处理需要的数据
		windowSize = round(1.0 * audio.sampleRate / goalRate);
		windowNum = ceil(1.0 * audio.sampleNum / windowSize);
		int fftSize = windowSize, processNum = 0, percentNow = 0;
		while ((fftSize & (-fftSize)) != fftSize) fftSize++;
		float* lastSpectrum = NULL; float* specturm = NULL;
		spData = new float[windowNum];
        rate = goalRate;

        // 绘制加载进度条
        if (processBar) {
            cout << "Converting to Specturm Data: ";
            cout << "["; for (int i = 0; i < partNum; i++) cout << "-"; cout << "] 0/" << partNum;
        }
        cout.flush();

        // 对每个窗口做 fft 运算
		for (int i = 0; i < audio.sampleNum; i += windowSize) {
			float* tmp = new float[fftSize];
			for (int j = 0; j < fftSize; j++) tmp[j] = 0;
			for (int j = 0; j < windowSize && i + j < audio.sampleNum; j++)
				tmp[j] = (audio.fetchSampleData(i + j).first + audio.fetchSampleData(i + j).second) / audio.numChannels; // 声道合并
			comp* complex = complexize(tmp, tmp + fftSize);
			fft(complex, complex + fftSize);
			specturm = floatize(complex, complex + fftSize);
			for (int j = 0; j < fftSize; j++) specturm[j] *= 2.0 / fftSize;
			float newSF = 0;
			for (int j = 0; j < fftSize; j++) {
				float value = specturm[j] - (i != 0 ? lastSpectrum[j] : 0);
				if (value >= 0) newSF += value; 
			}
			free(tmp); free(complex); free(lastSpectrum);
			lastSpectrum = specturm;
			spData[processNum] = newSF; processNum++;
			if (processBar && processNum * partNum / windowNum != percentNow) {
                // 绘制进度条
                int siz = partNum - percentNow + 3 + to_string(percentNow).size() + to_string(partNum).size();
                for (int j = 0; j < siz; j++) cout << "\b";
				percentNow = processNum * partNum / windowNum;

                cout << "#";
                for (int j = 0; j < partNum - percentNow; j++) cout << "-";
                cout << "] " << percentNow << "/" << partNum;
                cout.flush();
			}
		}

        if (processBar) cout << endl;
        cout.flush();
        _full = true;
    };

    void drawGraph(const char* dest, float multiRate = 100) {
        if (!_full) {
            cerr << "Failed to draw specturm graph: no specturm data was filled" << endl;
            abort();
        }

        ofstream fout(dest);
        if (!fout) {
            cerr << "Failed to draw specturm graph: couldn't open file" << endl;
            abort();
        }

        // 输出图像
		fout << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">" << endl;
        for (int i = 1; i < windowNum; i++) {
            fout << "<line x1=\"" << i - 1 << "\" y1=\"" << spData[i - 1] * multiRate
                 << "\" x2=\"" << i << "\" y2=\"" << spData[i] * multiRate
                 << "\" style=\"stroke: rgb(255, 0, 0);\"/>" << endl;
        }
		fout << "</svg>" << endl;
        fout.close();
    }

    bool isFull() {
        return _full;
    }

    float fetchSpecturmData(int dataId) {
        if (dataId < 0 || dataId >= windowNum) {
            cerr << "Failed to fetch specturm data: out of range" << endl;
            abort();
        }
        return spData[dataId];
    }
    float* fetchSpecturmData() {
        return spData;
    }
};

class beatData {
    private:
    bool _full = false;
    bool _calc = false;

    public:
    int beatNum = 0;
    float* beatTime;
    int rate = 0;
    int averageSize = 10;
    float multiplier = 1.5;
    float filterLimit = 0.2;
    float* spData;
    float* threshold;
    float* prunnedSpecturm;
    int windowNum;
    beatData(){};
    beatData(specturmData sp, int averageSize = 100, float multiplier = 1.5) {
        if (!sp.isFull()) {
            cerr << "Failed to calculate beat data: no specturm data was filled" << endl;
            abort();
        }
        windowNum = sp.windowNum;
        rate = sp.rate;
        spData = sp.fetchSpecturmData();
        this -> multiplier = multiplier;

        // 计算平均值
        threshold = new float[windowNum];
        for (int i = 0; i < windowNum; i++) {
            int beg = max(0ll, i - averageSize / 2);
            int end = min(windowNum - 1, i + averageSize / 2);
            float summarySample = 0;
            for (int j = beg; j <= end; j++) summarySample += spData[j];
            threshold[i] = summarySample / (end - beg + 1) * multiplier;
        }

        // 差值处理
        prunnedSpecturm = new float[windowNum];
		for (int i = 0; i < windowNum; i++) {
			if (threshold[i] <= spData[i]) prunnedSpecturm[i] = spData[i] - threshold[i];
			else prunnedSpecturm[i] = 0;
		}

        _full = true;
    }

    void onsetDetection(float filterLimit = 0.2) {
        if (!_full) {
            cerr << "Failed to execute onset detection algorithm: no specturm data was initialized" << endl;
            abort();
        }

        this -> filterLimit = filterLimit;
        
        // 初步检测节拍点个数
		float maxSpecturm = 0, whereSpecturm;
		for (int i = 0; i < windowNum; i++) {
			if (prunnedSpecturm[i]) {
				if (prunnedSpecturm[i] > maxSpecturm) whereSpecturm = i;
				maxSpecturm = max(maxSpecturm, prunnedSpecturm[i]);
			} else {
				if (maxSpecturm > filterLimit) beatNum++;
				maxSpecturm = 0;
			}
		}

        // 获得节拍点
		int solvedNum = 0;
		beatTime = new float[beatNum];
		maxSpecturm = 0;
		for (int i = 0; i < windowNum; i++) {
			if (prunnedSpecturm[i]) {
				if (prunnedSpecturm[i] > maxSpecturm) whereSpecturm = i;
				maxSpecturm = max(maxSpecturm, prunnedSpecturm[i]);
			} else {
				if (maxSpecturm > filterLimit) beatTime[solvedNum++] = 1.0 * whereSpecturm / rate;
				maxSpecturm = 0;
			}
		}

        _calc = true;
    }

    void drawGraph(const char* dest, float multiRate = 100) {
        if (!_full) {
            cerr << "Failed to draw threshold graph: no specturm data was initialized" << endl;
            abort();
        }

        ofstream fout(dest);
        if (!fout) {
            cerr << "Failed to draw threshold graph: couldn't open file" << endl;
            abort();
        }

        // 输出图像
		fout << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">" << endl;
        for (int i = 1; i < windowNum; i++) {
            fout << "<line x1=\"" << i - 1 << "\" y1=\"" << spData[i - 1] * multiRate
                 << "\" x2=\"" << i << "\" y2=\"" << spData[i] * multiRate
                 << "\" style=\"stroke: rgb(255, 0, 0);\"/>" << endl;
        }
        for (int i = 1; i < windowNum; i++) {
            fout << "<line x1=\"" << i - 1 << "\" y1=\"" << threshold[i - 1] * multiRate
                 << "\" x2=\"" << i << "\" y2=\"" << threshold[i] * multiRate
                 << "\" style=\"stroke: rgb(0, 255, 0);\"/>" << endl;
        }
        for (int i = 1; i < windowNum; i++) {
            fout << "<line x1=\"" << i - 1 << "\" y1=\"" << prunnedSpecturm[i - 1] * multiRate
                 << "\" x2=\"" << i << "\" y2=\"" << prunnedSpecturm[i] * multiRate
                 << "\" style=\"stroke: rgb(0, 0, 255);\"/>" << endl;
        }
		fout << "</svg>" << endl;
        fout.close();
    }

    float* fetchBeatTime() {
        if (!_calc) {
            cerr << "Failed to fetch beat data: onset detection wasn't executed" << endl;
            abort();
        }
        return beatTime;
    }

    int fetchBeatNum() {
        if (!_calc) {
            cerr << "Failed to fetch beat number: onset detection wasn't executed" << endl;
            abort();
        }
        return beatNum;
    }
};

#undef int
#endif
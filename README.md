# beat-detection-cplusplus

由 C++ 实现的 Beat Detection 算法，可用于识别歌曲内的节拍点。

目前算法版本为 v1.0，处于测试阶段，仅可识别 wav 格式的歌曲。通过测试的 wav 格式有: `44100Hz 16bit`，不保证其他格式能够正常运行。

参考代码 (即 `main.cpp`) 可以实现识别歌曲节拍点并输出为 BanG Dream! 格式的谱面文件，目前版本仅支持输出蓝键。

编译使用 GNU G++，不需要任何额外参数及依赖文件。

具体使用：

`wavfile::open(const char*)`: 打开一个 wav 文件，括号内参数为文件路径。

`wavfile::waveChannel()`: 将 wav 文件内的数据转化为浮点数据。

`wavfile::getSpectrum(int)`: 以传入的参数为目标频率，使用 STFT 算法将歌曲时域转化为频域。

`wavfile::beatDetection()`: 根据频域数据，识别歌曲的节拍点。

`wavfile::beatSize`: 识别到的节拍点个数。

`wavfile::beats`: 识别到的节拍点的时间信息，该信息除以传入的目标频率极为该节拍点的出现时间。

`wavfile::averageSize`: 差分处理时要取平均值的数据范围。

`wavfile::multiplier`: 差分处理时平均值扩增倍数，数值越大最后得到的节拍点越少。

`wavfile::limit`: 识别节拍点时的识别限制，大于该限制的节拍点才会被识别。

经测试，当 multiplier 为 1 时，一首三分钟的歌曲能识别出 2800 个节拍点，而为 1.5 时，仅能识别出 700 个节拍点。

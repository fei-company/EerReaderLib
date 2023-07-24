#include <chrono>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <fstream>
#include <tiffio.h>

#include "../ElectronCountedFramesDecompressor.h"

template <typename T>
void SetTag(std::shared_ptr<TIFF> tiffFile, uint32_t tag, T value)
{
    TIFFSetField(tiffFile.get(), tag, value);
}

std::shared_ptr<TIFF> OpenFile(const std::wstring &filepath) {
    auto tiffFile = TIFFOpenW(filepath.c_str(), "w");
    if (tiffFile == nullptr) {
        throw std::runtime_error("Could not write file");
    }
    return std::shared_ptr<TIFF>(tiffFile, [](TIFF *f) { TIFFClose(f); });
}

void SaveTiff(const std::wstring &filepath, uint32_t width, uint32_t height,
    std::vector<uint8_t> decompressedImage) {
    auto tiffFile = OpenFile(filepath);

    SetTag<uint32_t>(tiffFile, TIFFTAG_IMAGEWIDTH, width);
    SetTag<uint32_t>(tiffFile, TIFFTAG_IMAGELENGTH, height);
    SetTag<uint32_t>(tiffFile, TIFFTAG_ROWSPERSTRIP, height);
    SetTag<uint16_t>(tiffFile, TIFFTAG_BITSPERSAMPLE, 8);
    SetTag<uint16_t>(tiffFile, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);

    SetTag<uint16_t>(tiffFile, TIFFTAG_SAMPLESPERPIXEL, 1);
    SetTag<float>(tiffFile, TIFFTAG_XRESOLUTION, 1);
    SetTag<float>(tiffFile, TIFFTAG_YRESOLUTION, 1);
    SetTag<uint16_t>(tiffFile, TIFFTAG_RESOLUTIONUNIT, RESUNIT_NONE);


    auto bufferSize = width * height;
    // http://www.libtiff.org/man/TIFFWriteEncodedStrip.3t.html
    TIFFWriteEncodedStrip(tiffFile.get(), 0, reinterpret_cast<void *>(decompressedImage.data()),
        bufferSize);

    // http://www.libtiff.org/man/TIFFWriteDirectory.3t.html
    TIFFWriteDirectory(tiffFile.get());
}


int main(int argc, char* argv[]) {
    std::string inputFile;
    std::wstring outputFile;
    int upscaleFactor;

    if (argc == 4) {
        inputFile = argv[1];
        outputFile = std::wstring(argv[2], argv[2] + strlen(argv[2]));
        upscaleFactor = std::stoi(argv[3]);
    }
    else if (argc == 3)
    {
        inputFile = argv[1];
        outputFile = std::wstring(argv[2], argv[2] + strlen(argv[2]));
        upscaleFactor = 1;
    }
    else if (argc == 2 && (strcmp(argv[1], "--help") == 0 ||
        strcmp(argv[1], "-h") == 0))
    {
        std::cout << "Provide two arguments:\n Argument 1: String with Full Path to the input eer file. \n Argument 2: String with the Full Path to the output tiff file. \n Argument 3: Upscale factor, a positive or negative Integer of factor 2 (Default value 1)" << std::endl;
        return EXIT_SUCCESS;
    }
    else {
        std::cout << "Not enough arguments provided. Use --help to understand the "
            "correct use.";
        return EXIT_FAILURE;
    }

    unsigned width;
    unsigned height;
    unsigned nrOfFrames;

    ElectronCountedFramesDecompressor decompressor(inputFile);
    decompressor.getSize(width, height, nrOfFrames);
    std::cout << "Width: " << width
        << " Height: " << height << std::endl;

    width *= upscaleFactor;
    height *= upscaleFactor;
    const auto bufferSize = width * height;
    std::vector<uint8_t> eerImage(bufferSize);

    std::cout << "Decoding EER Data ..." << std::endl;
    // First frame
    decompressor.decompressImage(eerImage.data(), upscaleFactor);
    // Rest of the frames
    for (unsigned int eerFrame = 1; eerFrame < nrOfFrames; eerFrame++) {
        decompressor.decompressImage_AddTo(eerImage.data(), upscaleFactor);
    }
    std::cout << "Number of electron hits: "
        << decompressor.getNElectronsCounted() << std::endl;
    std::cout << "Saving output as Tiff image ..." << std::endl;
    SaveTiff(outputFile + L".Tiff", width, height, eerImage);

    return EXIT_SUCCESS;
}
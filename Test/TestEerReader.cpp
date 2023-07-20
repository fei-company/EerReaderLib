#include <chrono>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <fstream>
#include <tiffio.h>
 
 #include "../ElectronCountedFramesDecompressor.h"

const uint16_t TIFFTAG_EER_ACQUISITION_METADATA = 65001;

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
              std::vector<uint16_t> decompressedImage) {
  auto tiffFile = OpenFile(filepath);

  SetTag<uint32_t>(tiffFile, TIFFTAG_IMAGEWIDTH, width);
  SetTag<uint32_t>(tiffFile, TIFFTAG_IMAGELENGTH, height);
  SetTag<uint32_t>(tiffFile, TIFFTAG_ROWSPERSTRIP, height);
  SetTag<uint16_t>(tiffFile, TIFFTAG_BITSPERSAMPLE, 16);
  SetTag<uint16_t>(tiffFile, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);

  SetTag<uint16_t>(tiffFile, TIFFTAG_SAMPLESPERPIXEL, 1);
  SetTag<float>(tiffFile, TIFFTAG_XRESOLUTION, 1);
  SetTag<float>(tiffFile, TIFFTAG_YRESOLUTION, 1);
  SetTag<uint16_t>(tiffFile, TIFFTAG_RESOLUTIONUNIT, RESUNIT_NONE);

  
  auto bufferSize = width * height;
  // http://www.libtiff.org/man/TIFFWriteEncodedStrip.3t.html
    TIFFWriteEncodedStrip(tiffFile.get(), 0, reinterpret_cast<void *>(decompressedImage.data()),
                                 bufferSize * 2);
    
  // http://www.libtiff.org/man/TIFFWriteDirectory.3t.html
    TIFFWriteDirectory(tiffFile.get());
}


 int main(int argc, char* argv[]) {
    std::string inputFile;
    std::wstring outputFile;
    if (argc == 3) 
    {
        inputFile = argv[1];
        outputFile = std::wstring(argv[2], argv[2] + strlen(argv[2]));
    } 
    else if (argc == 2 && (strcmp(argv[1], "--help") == 0 ||
                             strcmp(argv[1], "-h") == 0))
    {
        std::cout << "Provide two arguments:\n Argument 1: String with Full Path to the input eer file. \n Argument 2: String with the Full Path to the output tiff file" << std::endl;
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


    const auto bufferSize = width * height;
    std::vector<uint16_t> decompressedImage(bufferSize);
    std::vector<uint8_t> eerImage(bufferSize);

    for (unsigned int eerFrame = 0 ; eerFrame < nrOfFrames ; eerFrame++)
    {
        decompressor.decompressImage((eerImage.data()));
        uint16_t* out = decompressedImage.data(); 
        uint8_t *in = eerImage.data();
        for (unsigned int i = 0; i < bufferSize; i++)
        {
            *out += static_cast<uint8_t>(*in);
            ++out;
            ++in;
        }
    }

    SaveTiff(outputFile + L".Tiff", width, height, decompressedImage);

    std::cout << "Number of electron hits: " << decompressor.getNElectronsCounted() << std::endl;

    return EXIT_SUCCESS;
}
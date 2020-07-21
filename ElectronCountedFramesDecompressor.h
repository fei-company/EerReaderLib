 // Copyright (c) 2019 by FEI Company

#pragma once


#include <fstream>
#include <vector>
#include "FeiBitStreamer.h"
#include "EerFile.h"


const unsigned nBitsRLE = 7; // 8 for old prototype EER files
const unsigned nSubPixBits = 2;
const bool no_bit_waste_on_overflow_code = true; // false for old prototype EER files.

const unsigned nBitsPerCode = nBitsRLE + 2*nSubPixBits;


const unsigned cameraSize = 4096; // NOTE: ONLY SQUARE NOW! and only FACLON
const unsigned superResolutionFactor = (1<<nSubPixBits); // constant now since it is hard coded in the FPGA compressor anyway

const unsigned totalSuperResolutionImSize = superResolutionFactor * cameraSize;

const unsigned gainImageSize = 4096;





struct ElectronPos
{
    uint16_t x;
    uint16_t y;
    
    ElectronPos(uint16_t x, uint16_t y) : x(x), y(y) {}
};


class ElectronCountedFramesDecompressor
{
public:

	// constructor for read mode (gets sizes and options from file header)
    ElectronCountedFramesDecompressor(const std::string& filename);

    void getSize(unsigned& x, unsigned& y, unsigned& z);
    unsigned getNFrames();
    size_t getFileSize() { return fsizeBytes; }
    unsigned getNElectronsCounted() { return nElectronsCounted; }

    ///read entire image of specified size.
	void decompressImage(uint8_t* p, int superFactor=1, int frameNumber = -1);
    void decompressImage_AddTo(uint8_t* p, int superFactor=1, int frameNumber = -1);

	void decompressImage(float* p, int superFactor=1, int frameNumber = -1);
    void decompressImage_AddTo(float* p, int superFactor=1, int frameNumber = -1);

    //unsigned nElectronFrameUpperLimit();
    unsigned nElectronFractionUpperLimit(int frameStart, int frameStop);
    unsigned decompressCoordinateList(ElectronPos* pList, int frameNumber = -1);

    //void getNormalizedSubPixHist(float* r);


    ~ElectronCountedFramesDecompressor();


private:
    typedef uint64_t BitStreamWordType;

    void prepareRead();
    BitStreamer prepareFrameRead(int frameNumber = -1);
    void finalizeFrameRead(BitStreamer& myBitStreamer);
    void prepareCodec();

    void createIndex(); // creates the index on-the-fly for a headerless file.

    std::fstream fh;    // used in ECC mode 
    std::unique_ptr<Fei::Acquisition::EerReader::EerFile> eerFile; // used in TIFF mode

	// for ecc
    std::vector<uint64_t> frameStartPointers;
	std::vector<BitStreamWordType> globalBuffer; // only used in headerless mode.
    
    // for tiff
    std::vector< std::vector<unsigned char> > frameBuffers; // only used in headerless mode.

    int _frameCounter;

	unsigned nx, ny, nFrames;
    size_t fsizeBytes;
    unsigned nElectronsCounted;

    bool tiffMode;

    float subPixCounts[16];

};

 // Copyright (c) 2019 by FEI Company

#pragma once


#include <fstream>
#include <vector>
#include "FeiBitStreamer.h"
#include "EerFile.h"

const bool g_no_bit_waste_on_overflow_code = true; // false for old prototype EER files.

struct ElectronPos
{
    uint16_t x;
    uint16_t y;
    
    ElectronPos(uint16_t x, uint16_t y) : x(x), y(y) {}
};

struct EerFrameSettings
{
    uint32_t width;
    uint32_t lenght;

    uint16_t rleBits;
    uint16_t horzSubBits;
    uint16_t vertSubBits;

    uint16_t bitsPerCode;
    uint16_t horizontalSubBitsOffset;
    uint16_t widthBitsOffset;

    explicit EerFrameSettings(const Fei::Acquisition::EerReader::EerFrame* eerFrame);

    // use for ECC mode
    EerFrameSettings();
};

class ElectronCountedFramesDecompressor
{
public:

	// constructor for read mode (gets sizes and options from file header)
    ElectronCountedFramesDecompressor(const std::string& filename);

    void getSize(unsigned& x, unsigned& y, unsigned& z);
    unsigned getNFrames();
    size_t getFileSize() { return m_fsizeBytes; }
    unsigned getNElectronsCounted() { return m_nElectronsCounted; }

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

    std::fstream m_fh;    // used in ECC mode
    std::unique_ptr<Fei::Acquisition::EerReader::EerFile> m_eerFile; // used in TIFF mode

	// for ecc
    std::vector<uint64_t> m_frameStartPointers;
	std::vector<BitStreamWordType> m_globalBuffer; // only used in headerless mode.
    
    // for tiff
    std::vector< std::vector<unsigned char> > m_frameBuffers; // only used in headerless mode.

    int m_frameCounter;

    EerFrameSettings m_eerFrameSettings;
    unsigned m_nFrames;
    size_t m_fsizeBytes;
    unsigned m_nElectronsCounted;

    bool m_tiffMode;

    float m_subPixCounts[16];

};

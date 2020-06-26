// Copyright (c) 2019 by FEI Company

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cassert>
#include <stdexcept>
#include <cstring>
#include <cstdlib>
#include "stddef.h"
#include "stdint.h"

#include "ElectronCountedFramesDecompressor.h"




/*****************************
 * internally used functions *
 *****************************/

// purely for debugging; WILL FAIL IN CASE OF MULTITHREAD
//#define SUBPIXHIST

#define COMPRESSIONMODE 6

const unsigned headerPrefixLen=20;
const char* headerPrefix="ThermoFisherECCompr"; // 19 char + \0

const unsigned DEFAULTBUFFERSIZE=4096*4096/16/8; // assuming worst case 16x compression

const char HWfooterStringOK[] = "ThermoFisherECComprOK000";
const char HWfooterStringErr[] = "ThermoFisherECComprERR00";
const unsigned int HWfooterLength = 24;

enum footerCheckResult { footerOK, footerError, footerInvalid };


footerCheckResult CheckFooter(BitStreamWordType* data)
{
    char* footChar = reinterpret_cast<char*>(data);
    int footOKCount = 0;
    int footErrCount = 0;
    for (int i=0; i<HWfooterLength;++i)
    {
        if (footChar[i]==HWfooterStringOK[i])
            footOKCount++;
        if (footChar[i]==HWfooterStringErr[i])
            footErrCount++;
    }
    return (footOKCount==HWfooterLength)? footerOK : ((footErrCount==HWfooterLength)? footerError : footerInvalid);
}


#ifdef SUBPIXHIST
unsigned TSThist[16];
#endif

template<int upBits>
struct electronSetterDelta
{
    static const unsigned widthIn = 4096;
    static const unsigned widthOut = (widthIn << nSubPixBits) >> (nSubPixBits-upBits);
    
    inline static void apply(uint8_t* p, unsigned posX, unsigned posY, unsigned subPixX, unsigned subPixY)
    {
        unsigned posXS = ((posX << nSubPixBits) + (subPixX^2)) >> ((int)nSubPixBits-upBits);
        unsigned posYS = ((posY << nSubPixBits) + (subPixY^2)) >> ((int)nSubPixBits-upBits); 
#ifdef SUBPIXHIST
        TSThist[(subPixX^2)+(subPixY^2)*4] ++;
#endif
        p[widthOut * posYS + posXS] = 1;

    }
};


template<int upBits>
struct electronAdderDelta
{
    static const unsigned widthIn = 4096;
    static const unsigned widthOut = (widthIn << nSubPixBits) >> (nSubPixBits-upBits);

    inline static void apply(uint8_t* p, unsigned posX, unsigned posY, unsigned subPixX, unsigned subPixY)
    {
        unsigned posXS = ((posX << nSubPixBits) + (subPixX^2)) >> ((int)nSubPixBits-upBits);
        unsigned posYS = ((posY << nSubPixBits) + (subPixY^2)) >> ((int)nSubPixBits-upBits); 
#ifdef SUBPIXHIST
        TSThist[(subPixX^2)+(subPixY^2)*4] ++;
#endif
        p[widthOut * posYS + posXS] ++;
    }
};





template <typename T>
inline void GetWeightsSpline3(T w, T weights[])
{
	weights[3] = (T)(1.0 / 6.0) * w * w * w;
	weights[0] = (T)(1.0 / 6.0) + (T)(1.0 / 2.0) * w * (w - (T)1.0) - weights[3];
	weights[2] = w + weights[0] - (T)2.0 * weights[3];
	weights[1] = (T)1.0 - weights[0] - weights[2] - weights[3];
}

static const unsigned splineSize = 4;
static float BSplineLUT[4][4];

static const float splineScaleFactor = 9.9867; // pi^2 ;-)

void createBSplineLUT()
{
    float splineF[4];
    for (int i=0; i<4;++i)
    {
        float w = ((float)(i))/4 + 0.125;
        //std::cout<<" W "<<w<<std::endl;
        GetWeightsSpline3(w, splineF);
        for (int j=0; j<4; ++j)
        {
            BSplineLUT[i][j] = splineF[j] * splineScaleFactor;
            //std::cout << ", "<<BSplineLUT[i][j] <<std::endl;
        }
    }
}


template<int upBits>
struct electronAdderBSpline
{
    static const unsigned widthIn = 4096;
    static const int widthOut = widthIn << upBits;
    
    inline static void apply(float* p, int posX, int posY, int subX, int subY)
    {
        int posXS = ((posX - ((subX>>1)&1)) << upBits) + (subX >> (nSubPixBits-upBits));
        int posYS = ((posY - ((subY>>1)&1)) << upBits) + (subY >> (nSubPixBits-upBits)); 
        
        subX = (subX << upBits) & 3; //repl MASK
        subY = (subY << upBits) & 3; //repl MASK
        int yStart = (posYS>0)? 0 : 1; // repl other splines sizes now not taken care of...
        int yStop = (posYS>widthOut-3)?  (widthOut+1-posYS): splineSize;
        for (int y=yStart; y<yStop; ++y)
        {
            int xStart = (posXS>0)? 0 : 1; // repl other splines sizes now not taken care of...
            int xStop = (posXS>widthOut-3)?  (widthOut+1-posXS): splineSize;
            int posCur = (posYS+y-1)*widthOut + posXS-1;
            for (int x=xStart; x<xStop; ++x)
            {
                p[posCur+x] += BSplineLUT[subY][y] * BSplineLUT[subX][x];
            }
        }
#ifdef SUBPIXHIST
        TSThist[(subX^2)+(subY^2)*4] ++;
#endif
    }
};





template<class operationFunc, typename ImType>
unsigned doDecompressImage(BitStreamer& myBitStreamer, ImType* p, unsigned w, unsigned h)
{
    static const int maxVal = ((1<<nBitsRLE)-1);
	//std::cout<<"TST"<<(superPosFunc::upBitsVal)<<", "<<w<<", "<<h<<std::endl;
    //std::cerr<<"superfactorr"<<operationFunc::widthOut<<std::endl;
	
    int N = (int)(w*h);

    int symbol = (int)myBitStreamer.getBits(nBitsPerCode);
    int value = symbol & maxVal;
    int outCount = value;
    unsigned nElect = 0;
	while (outCount<N)
	{
		if (value < maxVal)
		{
            int subPix = symbol>>nBitsRLE;
            operationFunc::apply(p, (outCount & 4095), (outCount >> 12), (subPix & 3), (subPix >> 2)); 
			++outCount;
            ++nElect;
		}
        else if (no_bit_waste_on_overflow_code) //this is constant expression so I assume compiler eliminates the check
        {
            myBitStreamer.rewind(2*nSubPixBits); // only for no-waste EER RLE implementation.
        }
		symbol = (int)myBitStreamer.getBits(nBitsPerCode);
		value = symbol & maxVal;
		outCount += value;
	}
    if (outCount != N)
    {
        std::ostringstream oss;
        oss << "ElectronCountedFramesDecompressor: invalid RLE decoding. resulting outCount = " << outCount << " while it should equal total # pixels = " << N;
        //throw std::range_error(oss.str());
        std::cerr << "Warning "<<oss.str()<<std::endl;
    }
    //std::cout << " nDecodedTST "<<nDecodedTST<<", nOvflTST "<<nOvflTST<<std::endl;
    return nElect;
}








/***********************************
 * ElectronCountedFramesDecompressor impl *
 ***********************************/


ElectronCountedFramesDecompressor::ElectronCountedFramesDecompressor(const std::string& filename)
	: nElectronsCounted(0)
{
    tiffMode = (filename.substr(filename.length()-4)!=".ecc");
    if (tiffMode)
    {
        std::cout<<"ElectronCountedFramesDecompressor: reading using TIFF-EER mode." << std::endl;
        eerFile.reset(new Fei::Acquisition::EerReader::EerFile(filename));
    }
    else
    {
        std::cout<<"ElectronCountedFramesDecompressor: reading using ECC mode." << std::endl;
        fh.open(filename.c_str(), std::ios_base::in | std::ios::binary);
        if (!fh.is_open())
        {
            printf("Error: input file cannot be opened!");
            return;
        }
    }
    prepareRead();
#ifdef SUBPIXHIST
    for (int i=0; i<16;++i) TSThist[i]=0;
#endif
    createBSplineLUT();
}


void ElectronCountedFramesDecompressor::getSize(unsigned& x, unsigned& y, unsigned& z)
{
    x = nx;
    y = ny;
    z = nFrames;
}

unsigned ElectronCountedFramesDecompressor::getNFrames()
{
    return nFrames;
}

void ElectronCountedFramesDecompressor::decompressImage(uint8_t* p, int superFactor, int frameNumber)
{
    BitStreamer myBitStreamer = prepareFrameRead(frameNumber);
    unsigned nxo = nx, nyo = ny;
    if (superFactor > 0) { nxo *= superFactor; nyo *= superFactor; }
    if (superFactor < 0) { nxo /= (-superFactor); nyo /= (-superFactor); }
	std::memset(p, 0, nxo*nyo*sizeof(uint8_t)); //optionally make it possible to skip this if you _know_ it is already OK.
    unsigned nElect = 0;
    switch(superFactor)
    {
        case -32: nElect = doDecompressImage<electronSetterDelta<-5> >(myBitStreamer, p, nx, ny); break;
        case -16: nElect = doDecompressImage<electronSetterDelta<-4> >(myBitStreamer, p, nx, ny); break;
        case -8: nElect = doDecompressImage<electronSetterDelta<-3> >(myBitStreamer, p, nx, ny); break;
        case -4: nElect = doDecompressImage<electronSetterDelta<-2> >(myBitStreamer, p, nx, ny); break;
        case -2: nElect = doDecompressImage<electronSetterDelta<-1> >(myBitStreamer, p, nx, ny); break;
        case 1: nElect = doDecompressImage<electronSetterDelta<0> >(myBitStreamer, p, nx, ny); break;
        case 2: nElect = doDecompressImage<electronSetterDelta<1> >(myBitStreamer, p, nx, ny); break;
        case 4: nElect = doDecompressImage<electronSetterDelta<2> >(myBitStreamer, p, nx, ny); break;
        default: 
            std::cerr<<"Super sampling factor must be 1, 2, 4, -2, -4, -8, -16, or -32!"<<std::endl;
            throw std::range_error("Super sampling factor must be 1, 2, 4, -2, -4, -8, -16, or -32!");
    }
    nElectronsCounted += nElect;
	finalizeFrameRead(myBitStreamer);
}



void ElectronCountedFramesDecompressor::decompressImage_AddTo(uint8_t* p, int superFactor, int frameNumber)
{
    //std::cerr<<"superfactorr"<<superFactor<<std::endl;
	BitStreamer myBitStreamer = prepareFrameRead(frameNumber);
    unsigned nElect = 0;
    switch(superFactor)
    {
        case -32: nElect = doDecompressImage<electronAdderDelta<-5> >(myBitStreamer, p, nx, ny); break;
        case -16: nElect = doDecompressImage<electronAdderDelta<-4> >(myBitStreamer, p, nx, ny); break;
        case -8: nElect = doDecompressImage<electronAdderDelta<-3> >(myBitStreamer, p, nx, ny); break;
        case -4: nElect = doDecompressImage<electronAdderDelta<-2> >(myBitStreamer, p, nx, ny); break;
        case -2: nElect = doDecompressImage<electronAdderDelta<-1> >(myBitStreamer, p, nx, ny); break;
        case 1: nElect = doDecompressImage<electronAdderDelta<0> >(myBitStreamer, p, nx, ny); break;
        case 2: nElect = doDecompressImage<electronAdderDelta<1> >(myBitStreamer, p, nx, ny); break;
        case 4: nElect = doDecompressImage<electronAdderDelta<2> >(myBitStreamer, p, nx, ny); break;
        default: throw std::range_error("Super sampling factor must be 1, 2, 4, -2, -4, -8, -16, or -32!");
    }
    nElectronsCounted += nElect;
	finalizeFrameRead(myBitStreamer);	
}


void ElectronCountedFramesDecompressor::decompressImage(float* p, int superFactor, int frameNumber)
{
	std::memset(p, 0, nx*ny*superFactor*superFactor*sizeof(float)); //optionally make it possible to skip this if you _know_ it is already OK.
    decompressImage_AddTo(p, superFactor, frameNumber);
}

void ElectronCountedFramesDecompressor::decompressImage_AddTo(float* p, int superFactor, int frameNumber)
{
	BitStreamer myBitStreamer = prepareFrameRead(frameNumber);
    unsigned nElect = 0;
    switch(superFactor)
    {
        case 1: nElect = doDecompressImage<electronAdderBSpline<0> >(myBitStreamer, p, nx, ny); break;
        case 2: nElect = doDecompressImage<electronAdderBSpline<1> >(myBitStreamer, p, nx, ny); break;
        case 4: nElect = doDecompressImage<electronAdderBSpline<2> >(myBitStreamer, p, nx, ny); break;
        default: throw std::range_error("Super sampling factor must be 1, 2, or 4 in B-Spline mode!");
    }
    // default: throw std::range_error("Super sampling factor must be 0 (BSpline), 1, 2, or 4!");
    nElectronsCounted += nElect;
    finalizeFrameRead(myBitStreamer);
}

/*unsigned ElectronCountedFramesDecompressor::nElectronFrameUpperLimit()
{
    unsigned nElectUpperLimitMax = 0;
    for (int i = 0; i < nFrames; ++i)
    {
        size_t frameSize = frameStartPointers[i + 1] - frameStartPointers[i];
        unsigned nElectUpperLimit = (frameSize * 8 + nBitsPerCode - 1) / nBitsPerCode;
        if (nElectUpperLimit > nElectUpperLimitMax)
            nElectUpperLimitMax = nElectUpperLimit;
    }
    return nElectUpperLimitMax; // this is really the upper limit. exactly correct only if no 255-symbols, and all 64 bit exactly filled
}*/

unsigned ElectronCountedFramesDecompressor::nElectronFractionUpperLimit(int frameStart, int frameStop)
{
    size_t fractionSize = 0;
    if (tiffMode)
    {
        for (int i=frameStart; i<frameStop;++i)
            fractionSize += frameBuffers[i].size();
    }
    else
    {
        fractionSize = frameStartPointers[frameStop] - frameStartPointers[frameStart];
    }
    unsigned nElectUpperLimit = (fractionSize * 8 + nBitsPerCode - 1) / nBitsPerCode;
    return nElectUpperLimit;
}

unsigned ElectronCountedFramesDecompressor::decompressCoordinateList(ElectronPos* pList, int frameNumber)
{
      static const int maxVal = ((1<<nBitsRLE)-1);

    frameNumber = frameNumber % nFrames; // nice for testing;
    BitStreamer myBitStreamer = prepareFrameRead(frameNumber);
    unsigned nElect = 0;
    int N = (int)(nx*ny);

    int symbol = (int)myBitStreamer.getBits(nBitsPerCode);
    int value = symbol & maxVal;
    int outCount = value;
    while (outCount<N)
    {
        if (value < maxVal)
        {
            int subPix = symbol >> nBitsRLE;
            pList->x = (((outCount & 4095) << nSubPixBits) | ((subPix & 3) ^ 2));
            pList->y = (((outCount >> 12) << nSubPixBits) | ((subPix >> 2) ^ 2));
            ++pList;
            ++outCount;
            ++nElect;

#ifdef SUBPIXHIST
			TSThist[((subPix & 3)^2)+((subPix >> 2)^2)*4]++;
#endif  
        }
        else  if (no_bit_waste_on_overflow_code)
        {
            myBitStreamer.rewind(2*nSubPixBits); // only for no-waste EER RLE implementation.
        }
        symbol = (int)myBitStreamer.getBits(nBitsPerCode);
        value = symbol & maxVal;
        outCount += value;
    }
    if (outCount != N)
    {
        std::ostringstream oss;
        oss << "ElectronCountedFramesDecompressor: invalid RLE decoding. resulting outCount = " << outCount << " while it should equal total # pixels + 1 = " << N + 1;
        //throw std::range_error(oss.str());
        std::cerr << "Warning " << oss.str() << std::endl;
    }
    //std::cout << " nDecodedTST "<<nDecodedTST<<", nOvflTST "<<nOvflTST<<std::endl;
    nElectronsCounted += nElect;
    finalizeFrameRead(myBitStreamer);	
    return nElect;
}


ElectronCountedFramesDecompressor::~ElectronCountedFramesDecompressor()
{
#ifdef SUBPIXHIST
    for (int i=0; i<4;++i) 
     /*   std::cout << "Sub-pix histX ["<<i<<"] = "<< TSThist[i] <<  std::endl;
    for (int i=0; i<4;++i) 
        std::cout << "Sub-pix histY ["<<i<<"] = "<< TSThist[i+4] <<  std::endl;*/
    std::cout << "TSThist = {\n";
    for (int i=0; i<4;++i)
    { 
        std::cout << "     ";
        for (int j=0; j<4;++j) 
            std::cout << TSThist[j+i*4]<<",\t";
        std::cout << "\n";
    }
    std::cout << "}" << std::endl;
#endif
}




void ElectronCountedFramesDecompressor::prepareRead()
{
    unsigned short compressionMode, submode;
    unsigned firstFramePtr;
    // read in string and check
    char headerId[headerPrefixLen];    
    if (tiffMode)
    {
        nx = 4096;
        ny = 4096;
        
        //iterate through all frames to count them
        
        nFrames = 0;
        fsizeBytes = 0;
        while (auto frame = eerFile->GetNextEerFrame())
        {
            nFrames++;
            frameBuffers.push_back(frame->GetEerData());
            fsizeBytes += frameBuffers.back().size();
        }
        std::cout << "ElectronCountedFramesDecompressor::prepareRead: found "<<nFrames<<" frames in EER-TIFF file." << std::endl;

    }
    else
    {
        // get file size
        fh.seekg(0, std::ios::end);      // Place the file pointer at the end of file
        fsizeBytes = fh.tellg();
        fh.seekg(0);

        //raw file without a header.
        //std::cout<<"No FCD header found. Assuming it is a raw dumped compressed file ... \n"<<std::endl;
        // NOTE: this is the normal use case; headers are never written out in the current proto
        nx = 4096;
        ny = 4096;
        nFrames = 1;

        uint64_t streamSize = (fsizeBytes+sizeof(BitStreamWordType)-1) / sizeof(BitStreamWordType);
        globalBuffer.resize(streamSize);
        //std::cout<<"ssize..."<<streamSize<<std::endl;
        fh.seekg(0);
        fh.read((char*)(globalBuffer.data()), streamSize*sizeof(BitStreamWordType));
        createIndex();
    }
	_frameCounter = 0;
}

BitStreamer ElectronCountedFramesDecompressor::prepareFrameRead(int frameNumber)
{
    if (frameNumber >= 0)
	{
		_frameCounter = frameNumber;
	}
	//std::cout<<"ElectronCountedFramesDecompressor::prepareFrameRead : frame "<<_frameCounter<<", size "<<frameBuffers[_frameCounter].size() << std::endl;
    if (tiffMode)
    {
        return BitStreamer(reinterpret_cast<BitStreamWordType*>(frameBuffers[_frameCounter].data()));
    }
    else
    {
        uint64_t wordPos = frameStartPointers[_frameCounter] / sizeof(BitStreamWordType);
        //std::cout<<"prepareFramRead headerless: starting from bufPos "<<wordPos<<" to read frame " <<_frameCounter<<std::endl;
        return BitStreamer(globalBuffer.data() + wordPos);
    }
}

void ElectronCountedFramesDecompressor::finalizeFrameRead(BitStreamer& myBitStreamer)
{
    int nwr = myBitStreamer.nWordsRead();
    //uint64_t nWordsShouldBe = myBitStreamer.buffer[nwr];
    //std::cout<<" nwordsRead in bytes "<<nwr*8<<std::endl;
    _frameCounter++;
}

void ElectronCountedFramesDecompressor::createIndex()
{
    std::vector<int64_t> frameSizesReversed;
    int64_t reverseCount = globalBuffer.size();
    int nFramesOK = 0;
    int nFramesErr = 0;
    const int footerWords = HWfooterLength / sizeof(BitStreamWordType);
    while (reverseCount>0)
    {
        footerCheckResult r = CheckFooter(globalBuffer.data() + reverseCount - footerWords);
        if (r != footerInvalid)
        {
            if (r == footerOK) nFramesOK++;
            else if (r == footerError) nFramesErr++;
            int64_t frameSize = globalBuffer[reverseCount - footerWords - 1] + footerWords + 1;
            frameSizesReversed.push_back(frameSize);
            reverseCount -= frameSize;
        }
        else
            throw std::range_error("ElectronCountedFramesDecompressor::createIndex: Found incorrect footer string!");
    }
    nFrames = frameSizesReversed.size();
    frameStartPointers.resize(nFrames+1);
    frameStartPointers[0] = 0;
    uint64_t ptr = 0;
    for (int i=1; i<=nFrames; ++i)
    {
        ptr += frameSizesReversed[frameSizesReversed.size()-i] * sizeof(BitStreamWordType);
        frameStartPointers[i] = ptr;
    }
    std::cout << "ElectronCountedFramesDecompressor::createIndex: found "<<nFrames<<" frames; #OK = "<<nFramesOK<<", #Error = "<<nFramesErr << std::endl;
}


/*
 * Copyright (C) 2019 Thermo Fisher Scientific. Do not distribute.
 * This software is provided by Thermo Fisher Scientific to
 * - The Hospital for Sick Children, Toronto, Canada, under confidentiality conditions 
 *   described in the agreed “Equipment Evaluation Agreement” (signed November 2018),
 * - Structura Biotechnology Inc., Toronto, Canada, under confidentiality conditions 
 *   described in the agreed “Mutual Non-Disclosure Agreement” (signed February 2019).
 */

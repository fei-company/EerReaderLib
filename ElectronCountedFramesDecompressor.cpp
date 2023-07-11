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
    static const unsigned widthOut = (widthIn << g_nSubPixBits) >> (g_nSubPixBits-upBits);
    
    inline static void apply(uint8_t* p, unsigned posX, unsigned posY, unsigned subPixX, unsigned subPixY)
    {
        unsigned posXS = ((posX << g_nSubPixBits) + (subPixX^2)) >> ((int)g_nSubPixBits-upBits);
        unsigned posYS = ((posY << g_nSubPixBits) + (subPixY^2)) >> ((int)g_nSubPixBits-upBits); 
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
    static const unsigned widthOut = (widthIn << g_nSubPixBits) >> (g_nSubPixBits-upBits);

    inline static void apply(uint8_t* p, unsigned posX, unsigned posY, unsigned subPixX, unsigned subPixY)
    {
        unsigned posXS = ((posX << g_nSubPixBits) + (subPixX^2)) >> ((int)g_nSubPixBits-upBits);
        unsigned posYS = ((posY << g_nSubPixBits) + (subPixY^2)) >> ((int)g_nSubPixBits-upBits); 
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
        int posXS = ((posX - ((subX>>1)&1)) << upBits) + (subX >> (g_nSubPixBits-upBits));
        int posYS = ((posY - ((subY>>1)&1)) << upBits) + (subY >> (g_nSubPixBits-upBits)); 
        
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
    static const int maxVal = ((1<<g_nBitsRLE)-1);
	//std::cout<<"TST"<<(superPosFunc::upBitsVal)<<", "<<w<<", "<<h<<std::endl;
    //std::cerr<<"superfactorr"<<operationFunc::widthOut<<std::endl;
	
    int N = (int)(w*h);

    int symbol = (int)myBitStreamer.getBits(g_nBitsPerCode);
    int value = symbol & maxVal;
    int outCount = value;
    unsigned nElect = 0;
	while (outCount<N)
	{
		if (value < maxVal)
		{
            int subPix = symbol>>g_nBitsRLE;
            operationFunc::apply(p, (outCount & 4095), (outCount >> 12), (subPix & 3), (subPix >> 2)); 
			++outCount;
            ++nElect;
		}
        else if (g_no_bit_waste_on_overflow_code) //this is constant expression so I assume compiler eliminates the check
        {
            myBitStreamer.rewind(2*g_nSubPixBits); // only for no-waste EER RLE implementation.
        }
		symbol = (int)myBitStreamer.getBits(g_nBitsPerCode);
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
	: m_nElectronsCounted(0)
{
    m_tiffMode = (filename.substr(filename.length()-4)!=".ecc");
    if (m_tiffMode)
    {
        std::cout<<"ElectronCountedFramesDecompressor: reading using TIFF-EER mode." << std::endl;
        m_eerFile.reset(new Fei::Acquisition::EerReader::EerFile(filename));
    }
    else
    {
        std::cout<<"ElectronCountedFramesDecompressor: reading using ECC mode." << std::endl;
        m_fh.open(filename.c_str(), std::ios_base::in | std::ios::binary);
        if (!m_fh.is_open())
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
    x = m_nx;
    y = m_ny;
    z = m_nFrames;
}

unsigned ElectronCountedFramesDecompressor::getNFrames()
{
    return m_nFrames;
}

void ElectronCountedFramesDecompressor::decompressImage(uint8_t* p, int superFactor, int frameNumber)
{
    BitStreamer myBitStreamer = prepareFrameRead(frameNumber);
    unsigned nxo = m_nx, nyo = m_ny;
    if (superFactor > 0) { nxo *= superFactor; nyo *= superFactor; }
    if (superFactor < 0) { nxo /= (-superFactor); nyo /= (-superFactor); }
	std::memset(p, 0, nxo*nyo*sizeof(uint8_t)); //optionally make it possible to skip this if you _know_ it is already OK.
    unsigned nElect = 0;
    switch(superFactor)
    {
        case -32: nElect = doDecompressImage<electronSetterDelta<-5> >(myBitStreamer, p, m_nx, m_ny); break;
        case -16: nElect = doDecompressImage<electronSetterDelta<-4> >(myBitStreamer, p, m_nx, m_ny); break;
        case -8: nElect = doDecompressImage<electronSetterDelta<-3> >(myBitStreamer, p, m_nx, m_ny); break;
        case -4: nElect = doDecompressImage<electronSetterDelta<-2> >(myBitStreamer, p, m_nx, m_ny); break;
        case -2: nElect = doDecompressImage<electronSetterDelta<-1> >(myBitStreamer, p, m_nx, m_ny); break;
        case 1: nElect = doDecompressImage<electronSetterDelta<0> >(myBitStreamer, p, m_nx, m_ny); break;
        case 2: nElect = doDecompressImage<electronSetterDelta<1> >(myBitStreamer, p, m_nx, m_ny); break;
        case 4: nElect = doDecompressImage<electronSetterDelta<2> >(myBitStreamer, p, m_nx, m_ny); break;
        default: 
            std::cerr<<"Super sampling factor must be 1, 2, 4, -2, -4, -8, -16, or -32!"<<std::endl;
            throw std::range_error("Super sampling factor must be 1, 2, 4, -2, -4, -8, -16, or -32!");
    }
    m_nElectronsCounted += nElect;
	finalizeFrameRead(myBitStreamer);
}



void ElectronCountedFramesDecompressor::decompressImage_AddTo(uint8_t* p, int superFactor, int frameNumber)
{
    //std::cerr<<"superfactorr"<<superFactor<<std::endl;
	BitStreamer myBitStreamer = prepareFrameRead(frameNumber);
    unsigned nElect = 0;
    switch(superFactor)
    {
        case -32: nElect = doDecompressImage<electronAdderDelta<-5> >(myBitStreamer, p, m_nx, m_ny); break;
        case -16: nElect = doDecompressImage<electronAdderDelta<-4> >(myBitStreamer, p, m_nx, m_ny); break;
        case -8: nElect = doDecompressImage<electronAdderDelta<-3> >(myBitStreamer, p, m_nx, m_ny); break;
        case -4: nElect = doDecompressImage<electronAdderDelta<-2> >(myBitStreamer, p, m_nx, m_ny); break;
        case -2: nElect = doDecompressImage<electronAdderDelta<-1> >(myBitStreamer, p, m_nx, m_ny); break;
        case 1: nElect = doDecompressImage<electronAdderDelta<0> >(myBitStreamer, p, m_nx, m_ny); break;
        case 2: nElect = doDecompressImage<electronAdderDelta<1> >(myBitStreamer, p, m_nx, m_ny); break;
        case 4: nElect = doDecompressImage<electronAdderDelta<2> >(myBitStreamer, p, m_nx, m_ny); break;
        default: throw std::range_error("Super sampling factor must be 1, 2, 4, -2, -4, -8, -16, or -32!");
    }
    m_nElectronsCounted += nElect;
	finalizeFrameRead(myBitStreamer);	
}


void ElectronCountedFramesDecompressor::decompressImage(float* p, int superFactor, int frameNumber)
{
	std::memset(p, 0, m_nx*m_ny*superFactor*superFactor*sizeof(float)); //optionally make it possible to skip this if you _know_ it is already OK.
    decompressImage_AddTo(p, superFactor, frameNumber);
}

void ElectronCountedFramesDecompressor::decompressImage_AddTo(float* p, int superFactor, int frameNumber)
{
	BitStreamer myBitStreamer = prepareFrameRead(frameNumber);
    unsigned nElect = 0;
    switch(superFactor)
    {
        case 1: nElect = doDecompressImage<electronAdderBSpline<0> >(myBitStreamer, p, m_nx, m_ny); break;
        case 2: nElect = doDecompressImage<electronAdderBSpline<1> >(myBitStreamer, p, m_nx, m_ny); break;
        case 4: nElect = doDecompressImage<electronAdderBSpline<2> >(myBitStreamer, p, m_nx, m_ny); break;
        default: throw std::range_error("Super sampling factor must be 1, 2, or 4 in B-Spline mode!");
    }
    // default: throw std::range_error("Super sampling factor must be 0 (BSpline), 1, 2, or 4!");
    m_nElectronsCounted += nElect;
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
    if (m_tiffMode)
    {
        for (int i=frameStart; i<frameStop;++i)
            fractionSize += m_frameBuffers[i].size();
    }
    else
    {
        fractionSize = m_frameStartPointers[frameStop] - m_frameStartPointers[frameStart];
    }
    unsigned nElectUpperLimit = (fractionSize * 8 + g_nBitsPerCode - 1) / g_nBitsPerCode;
    return nElectUpperLimit;
}

unsigned ElectronCountedFramesDecompressor::decompressCoordinateList(ElectronPos* pList, int frameNumber)
{
      static const int maxVal = ((1<<g_nBitsRLE)-1);

    frameNumber = frameNumber % m_nFrames; // nice for testing;
    BitStreamer myBitStreamer = prepareFrameRead(frameNumber);
    unsigned nElect = 0;
    int N = (int)(m_nx*m_ny);

    int symbol = (int)myBitStreamer.getBits(g_nBitsPerCode);
    int value = symbol & maxVal;
    int outCount = value;
    while (outCount<N)
    {
        if (value < maxVal)
        {
            int subPix = symbol >> g_nBitsRLE;
            pList->x = (((outCount & 4095) << g_nSubPixBits) | ((subPix & 3) ^ 2));
            pList->y = (((outCount >> 12) << g_nSubPixBits) | ((subPix >> 2) ^ 2));
            ++pList;
            ++outCount;
            ++nElect;

#ifdef SUBPIXHIST
			TSThist[((subPix & 3)^2)+((subPix >> 2)^2)*4]++;
#endif  
        }
        else  if (g_no_bit_waste_on_overflow_code)
        {
            myBitStreamer.rewind(2*g_nSubPixBits); // only for no-waste EER RLE implementation.
        }
        symbol = (int)myBitStreamer.getBits(g_nBitsPerCode);
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
    m_nElectronsCounted += nElect;
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
    if (m_tiffMode)
    {
        m_nx = 4096;
        m_ny = 4096;
        
        //iterate through all frames to count them
        
        m_nFrames = 0;
        m_fsizeBytes = 0;
        while (auto frame = m_eerFile->GetNextEerFrame())
        {
            m_nFrames++;
            m_frameBuffers.push_back(frame->GetEerData());
            m_fsizeBytes += m_frameBuffers.back().size();
        }
        std::cout << "ElectronCountedFramesDecompressor::prepareRead: found "<<m_nFrames<<" frames in EER-TIFF file." << std::endl;

    }
    else
    {
        // get file size
        m_fh.seekg(0, std::ios::end);      // Place the file pointer at the end of file
        m_fsizeBytes = m_fh.tellg();
        m_fh.seekg(0);

        //raw file without a header.
        //std::cout<<"No FCD header found. Assuming it is a raw dumped compressed file ... \n"<<std::endl;
        // NOTE: this is the normal use case; headers are never written out in the current proto
        m_nx = 4096;
        m_ny = 4096;
        m_nFrames = 1;

        uint64_t streamSize = (m_fsizeBytes+sizeof(BitStreamWordType)-1) / sizeof(BitStreamWordType);
        m_globalBuffer.resize(streamSize);
        //std::cout<<"ssize..."<<streamSize<<std::endl;
        m_fh.seekg(0);
        m_fh.read((char*)(m_globalBuffer.data()), streamSize*sizeof(BitStreamWordType));
        createIndex();
    }
	m_frameCounter = 0;
}

BitStreamer ElectronCountedFramesDecompressor::prepareFrameRead(int frameNumber)
{
    if (frameNumber >= 0)
	{
		m_frameCounter = frameNumber;
	}
	//std::cout<<"ElectronCountedFramesDecompressor::prepareFrameRead : frame "<<_frameCounter<<", size "<<frameBuffers[_frameCounter].size() << std::endl;
    if (m_tiffMode)
    {
        return BitStreamer(reinterpret_cast<BitStreamWordType*>(m_frameBuffers[m_frameCounter].data()));
    }
    else
    {
        uint64_t wordPos = m_frameStartPointers[m_frameCounter] / sizeof(BitStreamWordType);
        //std::cout<<"prepareFramRead headerless: starting from bufPos "<<wordPos<<" to read frame " <<_frameCounter<<std::endl;
        return BitStreamer(m_globalBuffer.data() + wordPos);
    }
}

void ElectronCountedFramesDecompressor::finalizeFrameRead(BitStreamer& myBitStreamer)
{
    int nwr = myBitStreamer.nWordsRead();
    //uint64_t nWordsShouldBe = myBitStreamer.buffer[nwr];
    //std::cout<<" nwordsRead in bytes "<<nwr*8<<std::endl;
    m_frameCounter++;
}

void ElectronCountedFramesDecompressor::createIndex()
{
    std::vector<int64_t> frameSizesReversed;
    int64_t reverseCount = m_globalBuffer.size();
    int nFramesOK = 0;
    int nFramesErr = 0;
    const int footerWords = HWfooterLength / sizeof(BitStreamWordType);
    while (reverseCount>0)
    {
        footerCheckResult r = CheckFooter(m_globalBuffer.data() + reverseCount - footerWords);
        if (r != footerInvalid)
        {
            if (r == footerOK) nFramesOK++;
            else if (r == footerError) nFramesErr++;
            int64_t frameSize = m_globalBuffer[reverseCount - footerWords - 1] + footerWords + 1;
            frameSizesReversed.push_back(frameSize);
            reverseCount -= frameSize;
        }
        else
            throw std::range_error("ElectronCountedFramesDecompressor::createIndex: Found incorrect footer string!");
    }
    m_nFrames = frameSizesReversed.size();
    m_frameStartPointers.resize(m_nFrames+1);
    m_frameStartPointers[0] = 0;
    uint64_t ptr = 0;
    for (int i=1; i<=m_nFrames; ++i)
    {
        ptr += frameSizesReversed[frameSizesReversed.size()-i] * sizeof(BitStreamWordType);
        m_frameStartPointers[i] = ptr;
    }
    std::cout << "ElectronCountedFramesDecompressor::createIndex: found "<<m_nFrames<<" frames; #OK = "<<nFramesOK<<", #Error = "<<nFramesErr << std::endl;
}

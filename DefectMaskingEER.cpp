/*
 * Copyright (C) 2019 Thermo Fisher Scientific. Do not distribute.
 */


#include "DefectMaskingEER.h"

const defectNeighborInfo_t bitMask_DefectNeighborTypeMask = 0xC0;
const defectNeighborInfo_t bitMask_Defect1DInterpol = 0x40;
const defectNeighborInfo_t bitMask_Defect2DInterpol = 0x80;
const defectNeighborInfo_t bitMask_DefectCorner = 0xC0;

const defectNeighborInfo_t bitMask_DefectUpDown = 0x20;
const defectNeighborInfo_t bitMask_DefectBack = 0x10;
const defectNeighborInfo_t bitMask_DefectRight = 0x00;
const defectNeighborInfo_t bitMask_DefectLeft = bitMask_DefectBack;
const defectNeighborInfo_t bitMask_DefectDown = bitMask_DefectUpDown;
const defectNeighborInfo_t bitMask_DefectUp = bitMask_DefectBack | bitMask_DefectUpDown;


const defectNeighborInfo_t bitMask_Distance = 0x0F;



const unsigned gainImageSizeFactor = totalSuperResolutionImSize/gainImageSize; // constant now since it is hard coded in the FPGA compressor anyway



// small convenience inline function. was only useful for unordered_maps impl. but that was was WAY slower so reverted
inline void _addDefectNeighborInfo(CameraDefectNeighborInformation& camDefectNeighborInfo, uint32_t idx, defectNeighborInfo_t spec)
{
    //camDefectNeighborInfo[idx] = DefectNeighborSpec(spec, gainImage[idx]); 
    camDefectNeighborInfo.neighborSpec[idx] = spec;
}

void CreateDefectNeighborInfoMap(const CameraDefects& def, const float* gainImage, CameraDefectNeighborInformation& camDefectNeighborInfo)
{
    camDefectNeighborInfo.gainImage = gainImage; // copy the pointer, simply/
	for (auto it = def.vLineDefects.begin(); it != def.vLineDefects.end(); ++it)
	{
		int xStart = it->begin;
		int xStop = it->end;
        
        defectNeighborInfo_t defLen = (defectNeighborInfo_t)(xStop - xStart);
        for (int y = 0; y<4096; ++y)
        { 
            for (int xx=xStart; xx<=xStop;++xx)
            {
                int idx = (y)*4096 + xx;
                _addDefectNeighborInfo(camDefectNeighborInfo, idx,  1); // 1 is used to mark a defect itself.
            } 
            int ofs = y*4096;
            if (xStart>0)
            {
                _addDefectNeighborInfo(camDefectNeighborInfo, 
                                       ofs + xStart - 1, bitMask_Defect1DInterpol | bitMask_DefectRight | defLen);
            }
            if (xStop<4095)
            {
                _addDefectNeighborInfo(camDefectNeighborInfo, 
                                       ofs + xStop + 1, bitMask_Defect1DInterpol | bitMask_DefectLeft | defLen);
            }
        }
	}
    // dont change order! first vlime, then Hline!
	for (auto it = def.hLineDefects.begin(); it != def.hLineDefects.end(); ++it)
	{
		int yStart = it->begin;
		int yStop = it->end;
        
        defectNeighborInfo_t defLen = (defectNeighborInfo_t)(yStop - yStart);
        for (int x = 0; x<4096; ++x)
        { 
            for (int yy=yStart; yy<=yStop;++yy)
            {
                int idx = (yy)*4096 + x;
                _addDefectNeighborInfo(camDefectNeighborInfo, idx,  1); // 1 is used to mark a defect itself.
            } 
            if (yStart>0)
            {
                int idx = (yStart-1)*4096 + x;
                defectNeighborInfo_t spec =
                         ((camDefectNeighborInfo.neighborSpec[idx] & bitMask_DefectNeighborTypeMask) && yStart>1)? 
                                bitMask_DefectCorner : bitMask_Defect1DInterpol;
                _addDefectNeighborInfo(camDefectNeighborInfo, 
                                       idx,  spec | bitMask_DefectDown | defLen);

            }
            if (yStop<4095)
            {
                int idx = (yStop+1)*4096 + x;
                defectNeighborInfo_t spec =
                         ((camDefectNeighborInfo.neighborSpec[idx] & bitMask_DefectNeighborTypeMask) && yStop<4094)? 
                                bitMask_DefectCorner : bitMask_Defect1DInterpol;
                _addDefectNeighborInfo(camDefectNeighborInfo, 
                                       idx,  spec | bitMask_DefectUp | defLen);
            }
        }
	}
	for (auto it = def.pixelDefects.begin(); it != def.pixelDefects.end(); ++it)
	{
		int x = it->x;
		int y = it->y;
        // std::cout << "X: "<< x << " - Y: " << y << std::endl;

        const defectNeighborInfo_t defLen = 0; //always for pixels. prep for areas.
        _addDefectNeighborInfo(camDefectNeighborInfo,  y*4096 + x,  1);
        if (y > 0)
            _addDefectNeighborInfo(camDefectNeighborInfo,  (y-1)*4096 + x,  bitMask_Defect2DInterpol | bitMask_DefectDown  | defLen);

        if (x > 0)
            _addDefectNeighborInfo(camDefectNeighborInfo,  y*4096 + x-1,    bitMask_Defect2DInterpol | bitMask_DefectRight | defLen);
        
        if (x < 4095)
            _addDefectNeighborInfo(camDefectNeighborInfo,  y*4096 + x+1,    bitMask_Defect2DInterpol | bitMask_DefectLeft  | defLen);

        if (y < 4095)
            _addDefectNeighborInfo(camDefectNeighborInfo,  (y+1)*4096 + x,  bitMask_Defect2DInterpol | bitMask_DefectUp    | defLen);
	}
	for (auto it = def.areaDefects.begin(); it != def.areaDefects.end(); ++it)
	{
    		defectNeighborInfo_t lx = (it->endX - it->beginX);
        defectNeighborInfo_t ly = (it->endY - it->beginY);
        
        for (int yy=it->beginY; yy<=it->endY;++yy)
            for (int xx=it->beginX; xx<=it->endX;++xx)
                _addDefectNeighborInfo(camDefectNeighborInfo,  yy*4096 + xx,  1);

        if ((it->beginY) > 0)
        {
            int y = it->beginY-1;
            for (int x=it->beginX; x<=it->endX; ++x)
                _addDefectNeighborInfo(camDefectNeighborInfo,  
                                       y*4096 + x, bitMask_Defect2DInterpol | bitMask_DefectDown | ly);
        }
        if ((it->endY) < 4095-1)
        {
            int y = it->endY+1;
            for (int x=it->beginX; x<=it->endX; ++x)
                _addDefectNeighborInfo(camDefectNeighborInfo,  
                                       y*4096 + x, bitMask_Defect2DInterpol | bitMask_DefectUp | ly);
        }
        if ((it->beginX) > 0)
        {
            int x = it->beginX-1;
            for (int y=it->beginY; y<=it->endY; ++y)
                _addDefectNeighborInfo(camDefectNeighborInfo,  
                                       y*4096 + x, bitMask_Defect2DInterpol | bitMask_DefectRight | lx);
        }
        if ((it->endX) < 4095-1)
        {
            int x = it->endX+1;
            for (int y=it->beginY; y<=it->endY; ++y)
                _addDefectNeighborInfo(camDefectNeighborInfo,  
                                       y*4096 + x,  bitMask_Defect2DInterpol | bitMask_DefectLeft | lx);
        }
    }
}


DefectElectronAdder::DefectElectronAdder() : 
    distHit(0.0, 1.0), distSubPix(0, 15)
{
    std::random_device rd;
    generator = std::default_random_engine(rd()); // for replacing hot-pixel value
}

unsigned DefectElectronAdder::execute(ElectronPos* pListPtr, unsigned nElect, const CameraDefectNeighborInformation& camDefectNeighborInfo)
{
    defectCounts.clear();
    // visit all electrons and add to the map if they are direct neighbors
    for (unsigned i=0; i<nElect;++i)
    {
        uint32_t eOfs = ((pListPtr[i].y >> nSubPixBits) << 12) | (pListPtr[i].x >> nSubPixBits);
        uint32_t eOfsGain = ((pListPtr[i].y / gainImageSizeFactor) * gainImageSize) | (pListPtr[i].x /gainImageSizeFactor);
        defectNeighborInfo_t di = camDefectNeighborInfo.neighborSpec[eOfs];
        if (di==1)
        {
           // electron at defect!! make it harmless
           pListPtr[i].x = 0xffff;
           pListPtr[i].y = 0xffff;
        }
        if (di & bitMask_DefectNeighborTypeMask)
        {
            uint16_t subPix = distSubPix(generator);
            pListPtr[i].x = (pListPtr[i].x & 0xfffc) | (subPix&3);
            pListPtr[i].y = (pListPtr[i].y & 0xfffc) | (subPix>>2);            
            float pcGain = camDefectNeighborInfo.gainImage[eOfsGain];
            //std::cout<<"pcGain"<<pcGain<<std::endl;
            defectNeighborInfo_t defDist = (di & bitMask_Distance)+1;
            int32_t defOfs = (di & bitMask_DefectUpDown)? 4096:1; 
            defOfs = (di&bitMask_DefectBack)? -defOfs : +defOfs;
            uint8_t mult = ((di & bitMask_DefectNeighborTypeMask) == bitMask_Defect2DInterpol)? 2 : 1;
            for (uint8_t k=1; k<=defDist;++k)
            {
                defectCounts[eOfs + k*defOfs] += pcGain * ((float)((defDist-k+1)))/(mult*(defDist+1)); // unordered_map creates zero-initialized element on-the-fly if not there yet. so convenient, yet confusing
            }
            if ((di & bitMask_DefectNeighborTypeMask) == bitMask_DefectCorner && eOfs >= defOfs && eOfs < 4096*4096+defOfs)
            {
                defectNeighborInfo_t diOrt = camDefectNeighborInfo.neighborSpec[eOfs - defOfs];
                // get the info of the horizontal defect. (by construction in CreateDefectNeighborInfoImage, the corner points now have the vertical one)
                // note; gain should not be obtained from the nieghbor!
                defectNeighborInfo_t defDistOrt = (diOrt & bitMask_Distance)+1;
                int32_t defOfsOrt = (diOrt & bitMask_DefectUpDown)? 4096:1; 
                defOfsOrt = (diOrt & bitMask_DefectBack)? -defOfsOrt : +defOfsOrt;
                for (uint8_t k=0; k<=defDist;++k)
                {
                    for (uint8_t m=1; m<=defDistOrt;++m)
                    {
                        defectCounts[eOfs + k*defOfs + m*defOfsOrt] += 
                            pcGain * ((float)((defDist-k+1)))/(defDist+1) * ((float)((defDistOrt-m+1)))/(defDistOrt+1);
                    }
                }
            }
        }
    }

    //  probalistic addition of electrons.
    unsigned plc = nElect;
    for (const auto &pair : defectCounts) 
    {
        uint32_t defectOfs = pair.first;
        float electProb = pair.second;
        if (distHit(generator) < electProb)
        {
            uint16_t subPix = distSubPix(generator);
            uint16_t posX = (((defectOfs & 4095) << nSubPixBits) + (subPix>>2));
            uint16_t posY = (((defectOfs >> 12) << nSubPixBits) + (subPix&3));
            //if ((defectOfs >> 12)!=2142)
            //    std::cerr<<"HUH "<<posY<<"; "<<(((defectOfs >> 12) << nSubPixBits) + (subPix>>2))<<" ; "<<subPix<<std::endl;
            pListPtr[plc++] = ElectronPos(posX, posY); // What with subpixel info? random?
        }
    }
    return plc;
}


SubpixelPositionRandomizer::SubpixelPositionRandomizer() : 
    distSubPix(0, 15)
{
    std::random_device rd;
    generator = std::default_random_engine(rd()); // for replacing hot-pixel value
}

void SubpixelPositionRandomizer::execute(ElectronPos* pListPtr, unsigned nElect)
{
    for (unsigned i=0; i<nElect;++i)
    {
        uint16_t subPix = distSubPix(generator);
        pListPtr[i].x = (pListPtr[i].x & 0xfffc) | (subPix>>2);  
        pListPtr[i].y = (pListPtr[i].y & 0xfffc) | (subPix&3);  
    }
}




/*
 * Copyright (C) 2019 Thermo Fisher Scientific. Do not distribute.
 * This software is provided by Thermo Fisher Scientific to
 * - The Hospital for Sick Children, Toronto, Canada, under confidentiality conditions 
 *   described in the agreed “Equipment Evaluation Agreement” (signed November 2018),
 * - Structura Biotechnology Inc., Toronto, Canada, under confidentiality conditions 
 *   described in the agreed “Mutual Non-Disclosure Agreement” (signed February 2019).
 */

#pragma once 

#include <vector>
#include <unordered_map>
#include <random>
#include "ElectronCountedFramesDecompressor.h"
#include "GainDefectCorrect.h"


typedef uint8_t defectNeighborInfo_t;

struct DefectNeighborSpec
{
	defectNeighborInfo_t defectInfo;
	float gainCoefficient;

    DefectNeighborSpec() {}
    DefectNeighborSpec(defectNeighborInfo_t t, float g): defectInfo(t), gainCoefficient(g) {}
};


struct CameraDefectNeighborInformation
{
    std::vector<defectNeighborInfo_t> neighborSpec;
    const float* gainImage; // not data owner, so be careful not to free the gainImage array too early.

    CameraDefectNeighborInformation(): neighborSpec(4096*4096,0) {}
};

void CreateDefectNeighborInfoMap(const CameraDefects& def, const float* gainImage, CameraDefectNeighborInformation& camDefectNeighborInfo);



struct DefectElectronAdder
{
    DefectElectronAdder();

    unsigned execute(ElectronPos* pListPtr, unsigned nElect, const CameraDefectNeighborInformation& camDefectNeighborInfo);


private:
    std::default_random_engine generator;
    std::uniform_real_distribution<float> distHit;
    std::uniform_int_distribution<uint16_t> distSubPix;

    std::unordered_map<uint32_t, float> defectCounts;
};


// FOR SUPERRESOLUTION TESTING PURPOSES
struct SubpixelPositionRandomizer
{
    SubpixelPositionRandomizer();

    unsigned execute(ElectronPos* pListPtr, unsigned nElec);


private:
    std::default_random_engine generator;
    std::uniform_int_distribution<uint16_t> distSubPix;
};




/*
 * Copyright (C) 2019 Thermo Fisher Scientific. Do not distribute.
 * This software is provided by Thermo Fisher Scientific to
 * - The Hospital for Sick Children, Toronto, Canada, under confidentiality conditions 
 *   described in the agreed “Equipment Evaluation Agreement” (signed November 2018),
 * - Structura Biotechnology Inc., Toronto, Canada, under confidentiality conditions 
 *   described in the agreed “Mutual Non-Disclosure Agreement” (signed February 2019).
 */

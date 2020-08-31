/*
 * Copyright (C) 2019 Thermo Fisher Scientific. Do not distribute.
 * This software is provided by Thermo Fisher Scientific to
 * - The Hospital for Sick Children, Toronto, Canada, under confidentiality conditions 
 *   described in the agreed “Equipment Evaluation Agreement” (signed November 2018),
 * - Structura Biotechnology Inc., Toronto, Canada, under confidentiality conditions 
 *   described in the agreed “Mutual Non-Disclosure Agreement” (signed February 2019).
 */

// gain correction for 4k image reconstruction  without any shift or doming model applied.
// automatically extracts defects from gain image acquired using ECC proto.
#pragma once

#include <vector>
#include <iostream>

struct DefectPixel
{
	DefectPixel(uint16_t x, uint16_t y) : x(x), y(y) {}

	uint16_t x;
	uint16_t y;
};

struct DefectLines
{
	DefectLines(uint16_t b, uint16_t e) : begin(b), end(e) {}

	uint16_t begin;
	uint16_t end;
};

struct DefectArea
{
	DefectArea(uint16_t bX, uint16_t bY, uint16_t eX, uint16_t eY) : beginX(bX), beginY(bY), endX(eX), endY(eY) {}

	uint16_t beginX, beginY;
	uint16_t endX, endY;
};



struct CameraDefects
{
	std::vector<DefectPixel> pixelDefects;
	std::vector<DefectLines> vLineDefects;
	std::vector<DefectLines> hLineDefects;
	std::vector<DefectArea> areaDefects;
    
    unsigned nPixelsDefect;
};

// note: modifies gainImage
CameraDefects PrepareGainImage_ExtractDefects(float* gainImage, unsigned gainImSize, unsigned w = 4096, unsigned h = 4096, double deadThreshold = 0.05);

template <typename T>
void DefectCorrect(T* img, const CameraDefects& def, unsigned w=4096, unsigned h=4096);


/*
 * Copyright (C) 2019 Thermo Fisher Scientific. Do not distribute.
 * This software is provided by Thermo Fisher Scientific to
 * - The Hospital for Sick Children, Toronto, Canada, under confidentiality conditions 
 *   described in the agreed “Equipment Evaluation Agreement” (signed November 2018),
 * - Structura Biotechnology Inc., Toronto, Canada, under confidentiality conditions 
 *   described in the agreed “Mutual Non-Disclosure Agreement” (signed February 2019).
 */

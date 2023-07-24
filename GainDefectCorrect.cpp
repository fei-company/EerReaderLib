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

#include <vector>
#include <iostream>
#include <fstream>
#include "GainDefectCorrect.h"


void reduce16k_to_4k(const float *gain16k,float *gain4k)
{
    for (int y=0; y<4096; ++y)
        for (int x=0; x<4096; ++x)
        {
            double r = 0;
            for(int sy=0; sy<4; ++sy)
                for (int sx=0; sx<4; ++sx)
                {
                    r += gain16k[(y*4+sy)*16384 + (x*4+sx)];
                }
            gain4k[y*4096+x] = static_cast<float>(r);
        }
    std::ofstream ofile("TST4Kred.bin", std::ios::binary);
	ofile.write((char*)(gain4k), sizeof(float)*4096*4096);
}


CameraDefects PrepareGainImage_ExtractDefects(float* gainImageFull, unsigned gainImSize, unsigned w , unsigned h , double deadThreshold)
{
    std::cerr<<"PrepareGainImage_ExtractDefects "<<gainImSize<<"; "<<w<<"--"<<h<<", #PIX "<<std::endl;

    static const float defectPixThreshold = 0.9f;
    static const float defectPixVisitedMark = -2.0f;
    static const float defectPixVisitedThresh = -1.0f; //must be > defectPixVisitedMark and <0 to do a safe float comparison

    int scaleFactor = 1;
    float* gainImage = gainImageFull;
    std::vector<float> gainImageReduced;
    std::cout<<"GAINSIZECHECK "<<gainImSize <<"<"<< w;
    if (gainImSize > w)
    {
    	scaleFactor =  gainImSize/w;
    	gainImageReduced.resize(w*h);
    	reduce16k_to_4k(gainImageFull, gainImageReduced.data());
    	gainImage = gainImageReduced.data();
    }

	CameraDefects r;
    r.nPixelsDefect = 0;
	double sum = 0;
	unsigned nZero = 0, nNonZero = 0;
	float *gptr = gainImage;
	for (unsigned y = 0; y < h; ++y)
	{
		for (unsigned x = 0; x < w; ++x, ++gptr)
		{
			if (!(*gptr <= defectPixThreshold && *gptr>defectPixVisitedThresh))
			{
				sum += *gptr;
				nNonZero++;
			}
			else
				nZero++;
		}
	}
	double mean = sum / nNonZero;

	double hotPixelThreshold = 400 * mean;

	std::cout << std::endl << "Hot pixel threshold: " << hotPixelThreshold << std::endl;
	double nRemovedPixels = 0;

	// extract dead pixels and create the gain image.
	auto deadThrAbs = deadThreshold * mean;
	gptr = gainImage;
	for (unsigned y = 0; y < h; ++y)
	{
		for (unsigned x = 0; x < w; ++x, ++gptr)
		{
			// hot pixel removal
            if (*gptr > hotPixelThreshold) {
            	*gptr = 0; 
            	nRemovedPixels ++;
            } 

			// pixel is dead if smaller than death pixel threshold and not a line defect.
			if (*gptr < defectPixThreshold && *gptr>defectPixVisitedThresh)
			{
                // determine the region size. Note: not fool proof for non-square sized defect areas, which are rare anyway
                unsigned xo = x+1;
                while (xo < w)
                {
                	//std::cout<<"HMM "<<x<<","<<y<<","<<xo<<" : "<< gainImage[y*w + xo]<<std::endl;
                    if (gainImage[y*w + xo] >= defectPixThreshold) break;
                    ++xo;
                }
                unsigned yo = y;
                bool yEndReached = false;
                while (!yEndReached)
                {
                    // mark the previous row.
                    for (unsigned i=x; i<xo; ++i)
                        gainImage[yo*w + i] = defectPixVisitedMark;
                    ++yo;
                    if  (yo >= h)
                        break;
                    // search the current row.
                    for (unsigned i=x; i<xo; ++i)
                        if (gainImage[yo*w + i] > defectPixThreshold) 
                        { yEndReached=true; break;}
                }
                std::cout<<"FOUND "<<x<<"--"<<xo<<" ; "<<y<<"--"<<yo<<", #PIX "<<((xo-x)*(yo-y))<<std::endl;
                r.nPixelsDefect += (xo-x)*(yo-y);
                // determine defect type
                if ((xo-x)==1 && (yo-y)==1)
                {
                    //if (x > 0 && y > 0 && x < w - 1 && y < h - 1) // for simplicity we forbid  dead pixels at border.
                    r.pixelDefects.push_back(DefectPixel(x, y));
                }
                else if ((xo-x) == w) // it is a horizontal line defect since it spans entire image width
                {
                    r.hLineDefects.push_back(DefectLines(y,yo-1)); 
                }
                else if ((yo-y) == h) // it is a horizontal line defect since it spans entire image width
                {
                    r.vLineDefects.push_back(DefectLines(x,xo-1));
                }
                else // it is an area defect
                {
                    r.areaDefects.push_back(DefectArea(x,y,xo-1,yo-1));
                }
			}
			/*else if (*gptr < deadThrAbs)
			{
				*gptr = 1.0; // intermediate case. just don't blow up this pixel.... pixel is semi-dead. TODO CHK
			}
			else
			{
				*gptr = mean / (*gptr);
			}*/
		}
	}

	std::cout << "Number of new dead pixels: " << nRemovedPixels << std::endl;


	gptr = gainImageFull;
	double p = mean / scaleFactor / scaleFactor;
	for (unsigned y = 0; y < gainImSize; ++y)
	{
		for (unsigned x = 0; x < gainImSize; ++x, ++gptr)
		{
			if ((*gptr)>0)
				*gptr = static_cast<float>(p) / (*gptr);
		}
	}


	std::cout << "PrepareGainImage_ExtractDefects: Total # electrons=" << sum << ", nPixelsOK=" << nNonZero << ", so total dose for gain image was " << mean << "e/pix. #death pixels=" << nZero << "\nDefects::" << std::endl;
	std::cout << "nPixelsDefect: " << r.nPixelsDefect << std::endl;
	std::cout << "#horizontal line defects: " << r.hLineDefects.size() << std::endl;
	std::cout << "#vertical line defects: " << r.vLineDefects.size() << std::endl;
	for (auto it = r.vLineDefects.begin(); it != r.vLineDefects.end(); ++it)
		std::cout << "defect vertical line at " << (it->begin) << " -- " << (it->end) << std::endl;
	for (auto it = r.hLineDefects.begin(); it != r.hLineDefects.end(); ++it)
		std::cout << "defect horizontal line at " << (it->begin) << " -- " << (it->end) << std::endl;
	std::cout << "#pixel defects: " << r.pixelDefects.size() << std::endl;
	std::cout << "#area defects: " << r.areaDefects.size() << std::endl;

    // finally, make all defects have gain 1.0
	for (auto it = r.vLineDefects.begin(); it != r.vLineDefects.end(); ++it)
		for (unsigned int y = 0; y<gainImSize; ++y)
			for (int x = scaleFactor*(it->begin); x<scaleFactor*((it->end)+1); ++x)
				gainImageFull[y*gainImSize + x] = 1.0;
	for (auto it = r.hLineDefects.begin(); it != r.hLineDefects.end(); ++it)
		for (int y = scaleFactor*(it->begin); y<scaleFactor*((it->end)+1); ++y)
			for (unsigned int x = 0; x<gainImSize; ++x)
					gainImageFull[y*gainImSize + x] = 1.0;
	for (auto it = r.pixelDefects.begin(); it != r.pixelDefects.end(); ++it)
		for (int i=0; i<scaleFactor;++i)
			for (int j=0; j<scaleFactor;++j)
				gainImageFull[((it->y)*scaleFactor+i)*gainImSize + ((it->x)*scaleFactor+j)] = 1.0;
	for (auto it = r.areaDefects.begin(); it != r.areaDefects.end(); ++it)
		for (int y = scaleFactor*it->beginY; y<scaleFactor*((it->endY)+1); ++y)
			for (int x = scaleFactor*it->beginX; x<=scaleFactor*((it->endX)+1); ++x)
				gainImageFull[y*gainImSize + x] = 1.0;
	return r;
}



template <typename T>
void DefectCorrect(T* img, const CameraDefects& def, unsigned w , unsigned h ) // now only for 4k x 4k
{
	for (auto it = def.vLineDefects.begin(); it != def.vLineDefects.end(); ++it)
	{
		auto xStart = it->begin;
		auto xStop = it->end;

		for (unsigned int y = 0; y<h; ++y)
		{
			T pixLeft = img[y*h + ((xStart>0) ? xStart - 1 : xStop + 1)];
			T pixRight = img[y*h + ((xStop<(w - 1)) ? xStop + 1 : xStart - 1)];
			for (int x = xStart; x <= xStop; ++x)
			{
				double coef = static_cast<double>(x - xStart + 1) / (xStop - xStart + 2);
				img[y*w + x] = static_cast<uint8_t>(pixLeft * (1 - coef) + pixRight * coef);
			}
		}
	}
	for (auto it = def.hLineDefects.begin(); it != def.hLineDefects.end(); ++it)
	{
		int yStart = it->begin;
		int yStop = it->end;

		for (unsigned int x = 0; x<w; ++x)
		{
			T pixLeft = img[x + ((yStart>0) ? yStart - 1 : yStop + 1) * 4096];
			T pixRight = img[x + ((yStop<4095) ? yStop + 1 : yStart - 1) * 4096];
			for (int y = yStart; y <= yStop; ++y)
			{
				double coef = (y - yStart + 1) / (yStop - yStart + 2);
				img[y*w + x] = static_cast<uint8_t>(pixLeft * (1 - coef) + pixRight * coef);
			}
		}
	}
	for (auto it = def.pixelDefects.begin(); it != def.pixelDefects.end(); ++it)
	{
		auto x = it->x;
		auto y = it->y;
		double sum = 0;
        int n = 4;
        if (y > 0) sum += img[(y-1)*w + x]; else --n;
        if (x > 0) sum += img[y*w + x - 1]; else --n;
        if (x < w-1) sum += img[y*w + x + 1]; else --n;
        if (y < h-1) sum += img[(y+1)*w + x]; else --n;
        img[y*w + x] = static_cast<uint8_t>(sum / n);
	}
	for (auto it = def.areaDefects.begin(); it != def.areaDefects.end(); ++it)
	{
		int xBef = ((it->beginX) > 0)? it->beginX - 1 : it->endX + 1;
		int xAft = ((it->beginX) < w-1)? it->endX + 1 : it->beginX - 1;
		int yBef = ((it->beginY) > 0)? it->beginY - 1 : it->endY + 1;
		int yAft = ((it->beginY) < h-1)? it->endY + 1 : it->beginY - 1;
        int nx = 2 * (it->endX - it->beginX + 2);
        int ny = 2 * (it->endY - it->beginY + 2);
        for (int y = it->beginY; y <= it->endY; ++y)
            for (int x = it->beginX; x <= it->endX; ++x)
            {
                double coefX = static_cast<double>(x - it->beginX + 1) / nx;
                double coefY = static_cast<double>(y - it->beginY + 1) / ny;
				img[y*w + x] = static_cast<uint8_t>((img[yBef*w + x] * (0.5 - coefY) + img[yAft*w + x] * coefY + img[y*w + xBef] * (0.5 - coefX) + img[y*w + xAft] * coefX));
            }
	}

}
// instantiations
template void DefectCorrect<float>(float* img, const CameraDefects& def, unsigned w, unsigned h);
template void DefectCorrect<uint8_t>(uint8_t* img, const CameraDefects& def, unsigned w, unsigned h);
template void DefectCorrect<uint16_t>(uint16_t* img, const CameraDefects& def, unsigned w, unsigned h);



/*
 * Copyright (C) 2019 Thermo Fisher Scientific. Do not distribute.
 * This software is provided by Thermo Fisher Scientific to
 * - The Hospital for Sick Children, Toronto, Canada, under confidentiality conditions 
 *   described in the agreed “Equipment Evaluation Agreement” (signed November 2018),
 * - Structura Biotechnology Inc., Toronto, Canada, under confidentiality conditions 
 *   described in the agreed “Mutual Non-Disclosure Agreement” (signed February 2019).
 */

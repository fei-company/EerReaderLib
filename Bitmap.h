// Copyright (c) 2019 by FEI Company
// All rights reserved. This file includes confidential and proprietary information of FEI Company.

#pragma once

#include <string>
#include <tiffio.h>
#include <vector>

namespace Fei {
namespace Acquisition {
namespace EerReader {


class Bitmap
{
public:
    Bitmap(TIFF* tiff);

    uint32_t GetWidth() const;
    uint32_t GetLength() const;

    std::vector<unsigned char> GetImageData() const;
    int GetBitsPerSample() const;

private:
    uint32_t m_imageWidth;
    uint32_t m_imageLength;
    std::vector<unsigned char> m_imageData;
    uint16_t m_bitsPerSample;
};

} //namespace EerReader
} //namespace Acquisition
} //namespace Fei

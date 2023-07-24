// Copyright (c) 2019 by FEI Company


#pragma once

#include <string>
#include <tiffio.h>
#include <vector>

namespace Fei {
namespace Acquisition {
namespace EerReader {


class EerFrame
{
public:
    EerFrame(TIFF* tiff);

    uint32_t GetWidth() const;
    uint32_t GetLength() const;
    uint16_t GetRleBits() const;
    uint16_t GetHorzSubBits() const;
    uint16_t GetVertSubBits() const;

    std::vector<unsigned char> GetEerData() const;
    int GetEncodingVersion() const;

private:
    uint32_t m_imageWidth;
    uint32_t m_imageLength;
    uint16_t m_rleBits;
    uint16_t m_horzSubBits;
    uint16_t m_vertSubBits;

    std::vector<unsigned char> m_eerData;
    int m_encodingVersion;
};

} //namespace EerReader
} //namespace Acquisition
} //namespace Fei

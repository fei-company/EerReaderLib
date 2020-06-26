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

    std::vector<unsigned char> GetEerData() const;
    int GetEncodingVersion() const;

private:
    uint32_t m_imageWidth;
    uint32_t m_imageLength;
    std::vector<unsigned char> m_eerData;
    int m_encodingVersion;
};

} //namespace EerReader
} //namespace Acquisition
} //namespace Fei

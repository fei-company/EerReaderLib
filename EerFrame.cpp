// Copyright (c) 2019 by FEI Company

#include "EerFrame.h"
#include <string>
#include <stdexcept>

namespace Fei {
namespace Acquisition {
namespace EerReader {

namespace {
    const uint16_t TIFF_COMPRESSION_EER_V0 = 65000;
    const uint16_t TIFF_COMPRESSION_EER_V1 = 65001;
}

EerFrame::EerFrame(TIFF * tiff)
{
    TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &m_imageWidth);
    TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &m_imageLength);
    auto compression = uint16_t(0);
    TIFFGetField(tiff, TIFFTAG_COMPRESSION, &compression);

    if (compression == TIFF_COMPRESSION_EER_V1)
        m_encodingVersion = 1;
    else if (compression == TIFF_COMPRESSION_EER_V0)
        m_encodingVersion = 0;
    else
        throw std::runtime_error("EerFrame: Unknown compression encoding " + std::to_string(compression));

    auto nrOfStrips = TIFFNumberOfStrips(tiff);
    for (uint32_t i = 0; i < nrOfStrips; i++)
    {
        auto previousSize = m_eerData.size();
        auto stripSize = TIFFRawStripSize(tiff, i);
        m_eerData.resize(previousSize + stripSize);

        TIFFReadRawStrip(tiff, i, m_eerData.data() + previousSize, stripSize);
    }
}

uint32_t EerFrame::GetWidth() const
{
    return m_imageWidth;
}

uint32_t EerFrame::GetLength() const
{
    return m_imageLength;
}

std::vector<unsigned char> EerFrame::GetEerData() const
{
    return m_eerData;
}

int EerFrame::GetEncodingVersion() const
{
    return m_encodingVersion;
}

} //namespace EerReader
} //namespace Acquisition
} //namespace Fei

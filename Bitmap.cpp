// Copyright (c) 2019 by FEI Company


#include "Bitmap.h"
#include <string>
#include <stdexcept>

namespace Fei {
namespace Acquisition {
namespace EerReader {


Bitmap::Bitmap(TIFF * tiff)
{
    TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &m_imageWidth);
    TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &m_imageLength);
    auto compression = uint16_t(0);
    TIFFGetField(tiff, TIFFTAG_COMPRESSION, &compression);
    TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &m_bitsPerSample);

    if (compression != COMPRESSION_NONE)
        throw std::runtime_error("Bitmap image: Unknown compression encoding " + std::to_string(compression));

    auto nrOfStrips = TIFFNumberOfStrips(tiff);
    for (uint32_t i = 0; i < nrOfStrips; i++)
    {
        auto previousSize = m_imageData.size();
        auto stripSize = TIFFRawStripSize(tiff, i);
        m_imageData.resize(previousSize + stripSize);

        TIFFReadRawStrip(tiff, i, m_imageData.data() + previousSize, stripSize);
    }
    auto imageSize = m_imageWidth * m_imageLength * (m_bitsPerSample >> 3);
    if (imageSize != m_imageData.size())
    {
        throw std::runtime_error("Bitmap image: image size on disk does not match image dimensions");
    }
}

uint32_t Bitmap::GetWidth() const
{
    return m_imageWidth;
}

uint32_t Bitmap::GetLength() const
{
    return m_imageLength;
}

std::vector<unsigned char> Bitmap::GetImageData() const
{
    return m_imageData;
}

int Bitmap::GetBitsPerSample() const
{
    return m_bitsPerSample;
}

} //namespace EerReader
} //namespace Acquisition
} //namespace Fei

// Copyright (c) 2019 by FEI Company
// All rights reserved. This file includes confidential and proprietary information of FEI Company.

#include "EerFile.h"
#include <iostream>
#include <stdexcept>


namespace Fei {
namespace Acquisition {
namespace EerReader {

namespace
{
    const uint16_t TIFF_COMPRESSION_EER_V0 = 65000;
    const uint16_t TIFF_COMPRESSION_EER_V1 = 65001;

    const uint16_t TIFFTAG_EER_ACQUISITION_METADATA = 65001;
    const uint16_t TIFFTAG_EER_FINAL_IMAGE_PROCESSING_METADATA = 65005;
    const uint16_t TIFFTAG_EER_FINAL_IMAGE_METADATA = 65006;

    void MyTIFFOutputError(const char* module, const char* fmt, va_list ap)
    {
        std::cout << module << ": ";
        vprintf(fmt, ap);
        std::cout << std::endl;
    }


    std::string GetFieldAsStringOrDefault(TIFF* tiff, ttag_t tag, const std::string& default)
    {
        char* data = nullptr;
        uint32_t count = 0;
        if (TIFFGetField(tiff, tag, &count, &data) != 1) return default;
        return std::string(data, count);
    }

    void TagExtender(TIFF *tif)
    {
        static const TIFFFieldInfo fieldInfo[] =
        {
            { TIFFTAG_EER_ACQUISITION_METADATA, TIFF_VARIABLE, TIFF_VARIABLE, TIFF_UNDEFINED, FIELD_CUSTOM, true, true, "Acquisition metadata" },
            { TIFFTAG_EER_FINAL_IMAGE_PROCESSING_METADATA, TIFF_VARIABLE, TIFF_VARIABLE, TIFF_UNDEFINED, FIELD_CUSTOM, true, true, "Final image processing metadata" },
            { TIFFTAG_EER_FINAL_IMAGE_METADATA, TIFF_VARIABLE, TIFF_VARIABLE, TIFF_UNDEFINED, FIELD_CUSTOM, true, true, "Final image metadata" },
        };

     TIFFMergeFieldInfo(tif, fieldInfo, sizeof(fieldInfo) / sizeof(fieldInfo[0]));
}
}

EerFile::EerFile(const std::string & filename)
{
    TIFFSetErrorHandler(MyTIFFOutputError);
    TIFFSetWarningHandler(MyTIFFOutputError);
    TIFFSetTagExtender(TagExtender);

    m_tiff = std::shared_ptr<TIFF>(TIFFOpen(filename.c_str(), "rl"), [](TIFF* tiffPtr){ if (tiffPtr) TIFFClose(tiffPtr); });
    if (!m_tiff)
        throw std::runtime_error("unable to open tiff file");

    if (!IsCurrentFrameEERCompressed())
    {
        m_finalImageBitmap = std::make_shared<Bitmap>(m_tiff.get());
    }
    else
    {
        m_finalImageBitmap = nullptr;
    }

    m_acquisitionMetadata = GetFieldAsStringOrDefault(m_tiff.get(), TIFFTAG_EER_ACQUISITION_METADATA, "");
    m_finalImageMetadata = GetFieldAsStringOrDefault(m_tiff.get(), TIFFTAG_EER_FINAL_IMAGE_METADATA, "");
}

EerFile::~EerFile()
{
}

std::unique_ptr<EerFrame> EerFile::GetNextEerFrame()
{
    std::unique_ptr<EerFrame> eerFrame;
    while (m_nextFrameAvailable && !eerFrame)
    {
        if (IsCurrentFrameEERCompressed())
            eerFrame = std::make_unique<EerFrame>(m_tiff.get());
        m_nextFrameAvailable = (TIFFReadDirectory(m_tiff.get()) == 1);
    }

    return std::move(eerFrame);
}

std::shared_ptr<Bitmap> EerFile::GetFinalImage()
{
    return m_finalImageBitmap;
}

std::string EerFile::GetAcquisitionMetadata() const
{
    return m_acquisitionMetadata;
}

std::string EerFile::GetFinalImageMetadata() const
{
    return m_finalImageMetadata;
}

bool EerFile::IsCurrentFrameEERCompressed()
{
    uint16_t compression;
    TIFFGetField(m_tiff.get(), TIFFTAG_COMPRESSION, &compression);
    return compression == TIFF_COMPRESSION_EER_V0 || compression == TIFF_COMPRESSION_EER_V1;
}

} //namespace EerReader
} //namespace Acquisition
} //namespace Fei

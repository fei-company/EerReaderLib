// Copyright (c) 2019 by FEI Company

#pragma once

#include <memory>
#include <string>
#include <tiffio.h>
#include "EerFrame.h"
#include "Bitmap.h"

namespace Fei {
namespace Acquisition {
namespace EerReader {

class EerFile
{
public:
    EerFile(const std::string& filename);
    ~EerFile();

    std::unique_ptr<EerFrame> GetNextEerFrame();
    std::shared_ptr<Bitmap> GetFinalImage();
    std::string GetAcquisitionMetadata() const;

private:
    bool IsCurrentFrameEERCompressed();

private:
    std::shared_ptr<TIFF> m_tiff;
    bool m_nextFrameAvailable = true;
    std::shared_ptr<Bitmap> m_finalImageBitmap;
};

} //namespace EerReader
} //namespace Acquisition
} //namespace Fei

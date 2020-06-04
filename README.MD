Version History
===============

  ------------------------------------------------------------------------------
  **Version**   **Date**      **Changes**                   **Author(s)**
  ------------- ------------- ----------------------------- --------------------
  1.0           30 Aug 2019   Initial version               P. Bootsma

  2.0           02 Dec 2019   Added 7-bit compression       A. Dumitrescu

  2.1           10 Dec 2019   Added Final image directory   M. Balasubramanian

  2.2           17 Jan 2020   Added acquisition metadata    A. Dumitrescu,
                                                            
                                                            P. Bootsma
  ------------------------------------------------------------------------------

References
==========

  **No**        **Title**                 **Author(s)**
  ------------- ------------------------- --------------------
  \[TIFF\]      TIFF Revision 6.0         Adobe Systems Inc.
  \[BIGTIFF\]   The BigTIFF File Format   Aware Systems

Introduction
============

Some cameras can produce EER encoded images. These images are stored in
a file format called the EER file format, identified by the extension
“.eer”. EER files are based on BigTIFF, but define a number of
extensions to this format to store all EER-related information in one
file.

The EER file writer, and also the exact file format are still in
development. This document describes the currently implemented file
format. The following features are planned for future inclusion in EER
files:

-   Gain reference data

-   Defect mask

-   Frame metadata

Note that these are the current plans, and that plans can change.

BigTIFF Container Format
========================

The EER file format is an extension of the BigTIFF file format defined
by Aware Systems \[BIGTIFF\]. BigTIFF is an extension of the regular
TIFF file format \[TIFF\] to support files of more than 4 gigabyte. EER
files add a custom compression scheme for the EER data and custom tags
to include the relevant metadata.

This document describes how EER data is stored within the BigTIFF
container format. Interpretation of the data is outside of the scope of
this document. The LibTIFF library \[http://simplesystems.org/libtiff\]
can read BigTIFF files, and therefore also EER files.

TIFF is a flexible container format by nature, and readers should adhere
to a set of requirements as specified in section 7 of the TIFF standard.
The following requirements are the most important ones to ensure future
extensibility of the file format:

-   Entries with unknown tag number or field type must be ignored.

-   Image file directories with an unexpected compression scheme must be
    skipped.

Extensions
==========

This section describes the extensions added to the BigTIFF container
format.

Compression schemes
-------------------

The EER file format adds one compression scheme to BigTIFF. This value
is used in the compression field (tag 259) to indicate that an IFD
contains EER compressed data.

  **Value**   **Comment**
  ----------- ---------------------
  65000       EER compressed data

Fields
------

The EER file format adds several fields to BigTIFF. These fields can be
present in IFDs as indicated later in this document.

  **Field**              **Tag**   **Type**    **Comment**
  ---------------------- --------- ----------- -------------------------------------------------------------------------------------
  Acquisition Metadata   65001     UNDEFINED   A non-zero-terminated string containing metadata concerning the entire acquisition.

Image File Directories
======================

An EER file contains one or more image file directories (IFDs). Every
IFD contains either an EER frame or other image data.

EER Frames
----------

An IFD containing an EER frame is identified by the compression scheme
of the image data:

  **Field**     **Tag**   **Type**   **Value**   **Comment**
  ------------- --------- ---------- ----------- ---------------------------
  Compression   259       SHORT      65000       EER compressed data 8-bit
  Compression   259       SHORT      65001       EER compressed data 7-bit

EER compressed data is stored like a bilevel image in a regular TIFF. It
is stored in one or more strips. If multiple strips are used, the entire
frame can be reconstructed by concatenating the data in the separate
strips.

An EER frame contains a number of required and optional fields. Unknown
or unexpected fields must be ignored.

### Required fields

  **Field**           **Tag**   **Type**               **Comment**
  ------------------- --------- ---------------------- ------------------------------------------------------------------------------------------------
  Image Width         256       SHORT or LONG          The number of columns in the frame
  Image Length        257       SHORT or LONG          The number of rows in the frame
  Compression         259       SHORT                  The compression scheme used; 65000/1 for EER data
  Rows Per Strip      278       SHORT or LONG          The number of rows in each strip (except possibly the last strip, which can contain less rows)
  Strip Offsets       273       SHORT, LONG or LONG8   For each strip, the byte offset of that strip
  Strip Byte Counts   279       SHORT, LONG or LONG8   For each strip, the number of bytes in that strip

### Optional fields

There are no optional fields specific to EER.

Final Image
-----------

An IFD containing the final pipeline image is included as part of the
EER file if a final image was recorded with the acquisition. When
included, the final image is always in the first IFD of the EER file.
The presence of a final image in the first IFD is identified by the
compression tag, which is set to no compression in case of a final
image, as opposed to EER frames. Currently, the EER file format only
supports EER frames and one optional final image frame.

### Required fields

  **Field**           **Tag**   **Type**               **Comment**
  ------------------- --------- ---------------------- ------------------------------------------------------------------------------------------------
  Image Width         256       SHORT or LONG          The number of columns in the frame
  Image Length        257       SHORT or LONG          The number of rows in the frame
  Compression         259       SHORT                  The compression scheme used; COMPRESSION\_NONE (1) for final image data
  Bits per sample     258       SHORT                  The number of bits per channel. For final image data, it is 16 bits per sample (2 bytes)
  Rows Per Strip      278       SHORT or LONG          The number of rows in each strip (except possibly the last strip, which can contain less rows)
  Strip Offsets       273       SHORT, LONG or LONG8   For each strip, the byte offset of that strip
  Strip Byte Counts   279       SHORT, LONG or LONG8   For each strip, the number of bytes in that strip

Metadata
========

In the first IFD metadata is written that applies to the entire
acquisition. The metadata is written as a tag with ID 65001. This tag
contains an XML string structured as follows:

&lt;metadata&gt;\
&lt;item name=”\[name\]” unit=”\[unit\]”&gt;\[value\]&lt;/item&gt;\
&lt;item name=”\[name2\]”&gt;\[value2\]&lt;/item&gt;\
…\
&lt;/metadata&gt;

The unit attribute may be omitted for values that do not have a unit.

The values present in the field are the following:

  **Name**            **Unit**    **Comments**
  ------------------- ----------- -------------------------------------------------------------------------
  numberOfFrames                  Total number of frames in the acquisition
  exposureTime        s           Exposure time of the acquisition
  totalDose           e/pixel     Average amount of electrons per pixel in the entire acquisition
  meanDoseRate        e/pixel/s   Average amount of electrons per pixel per second during the acquisition
  sensorImageWidth    pixels      Width of the individual images
  sensorImageHeight   pixels      Height of the individual images



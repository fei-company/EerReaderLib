Version History
===============

<table>
<thead>
<tr class="header">
<th><strong>Version</strong></th>
<th><strong>Date</strong></th>
<th><strong>Changes</strong></th>
<th><strong>Author(s)</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>1.0</td>
<td>30 Aug 2019</td>
<td>Initial version</td>
<td>P. Bootsma</td>
</tr>
<tr class="even">
<td>2.0</td>
<td>02 Dec 2019</td>
<td>Added 7-bit compression</td>
<td>A. Dumitrescu</td>
</tr>
<tr class="odd">
<td>2.1</td>
<td>10 Dec 2019</td>
<td>Added Final image directory</td>
<td>M. Balasubramanian</td>
</tr>
<tr class="even">
<td>2.2</td>
<td>17 Jan 2020</td>
<td>Added acquisition metadata</td>
<td><p>A. Dumitrescu,</p>
<p>P. Bootsma</p></td>
</tr>
</tbody>
</table>

References
==========

<table>
<thead>
<tr class="header">
<th><strong>No</strong></th>
<th><strong>Title</strong></th>
<th><strong>Author(s)</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>[TIFF]</td>
<td>TIFF Revision 6.0</td>
<td>Adobe Systems Inc.</td>
</tr>
<tr class="even">
<td>[BIGTIFF]</td>
<td>The BigTIFF File Format</td>
<td>Aware Systems</td>
</tr>
</tbody>
</table>

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

<table>
<thead>
<tr class="header">
<th><strong>Value</strong></th>
<th><strong>Comment</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>65000</td>
<td>EER compressed data</td>
</tr>
</tbody>
</table>

Fields
------

The EER file format adds several fields to BigTIFF. These fields can be
present in IFDs as indicated later in this document.

<table>
<thead>
<tr class="header">
<th><strong>Field</strong></th>
<th><strong>Tag</strong></th>
<th><strong>Type</strong></th>
<th><strong>Comment</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Acquisition Metadata</td>
<td>65001</td>
<td>UNDEFINED</td>
<td>A non-zero-terminated string containing metadata concerning the entire acquisition.</td>
</tr>
</tbody>
</table>

Image File Directories
======================

An EER file contains one or more image file directories (IFDs). Every
IFD contains either an EER frame or other image data.

EER Frames
----------

An IFD containing an EER frame is identified by the compression scheme
of the image data:

<table>
<thead>
<tr class="header">
<th><strong>Field</strong></th>
<th><strong>Tag</strong></th>
<th><strong>Type</strong></th>
<th><strong>Value</strong></th>
<th><strong>Comment</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Compression</td>
<td>259</td>
<td>SHORT</td>
<td>65000</td>
<td>EER compressed data 8-bit</td>
</tr>
<tr class="even">
<td>Compression</td>
<td>259</td>
<td>SHORT</td>
<td>65001</td>
<td>EER compressed data 7-bit</td>
</tr>
</tbody>
</table>

EER compressed data is stored like a bilevel image in a regular TIFF. It
is stored in one or more strips. If multiple strips are used, the entire
frame can be reconstructed by concatenating the data in the separate
strips.

An EER frame contains a number of required and optional fields. Unknown
or unexpected fields must be ignored.

### Required fields

<table>
<thead>
<tr class="header">
<th><strong>Field</strong></th>
<th><strong>Tag</strong></th>
<th><strong>Type</strong></th>
<th><strong>Comment</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Image Width</td>
<td>256</td>
<td>SHORT or LONG</td>
<td>The number of columns in the frame</td>
</tr>
<tr class="even">
<td>Image Length</td>
<td>257</td>
<td>SHORT or LONG</td>
<td>The number of rows in the frame</td>
</tr>
<tr class="odd">
<td>Compression</td>
<td>259</td>
<td>SHORT</td>
<td>The compression scheme used; 65000/1 for EER data</td>
</tr>
<tr class="even">
<td>Rows Per Strip</td>
<td>278</td>
<td>SHORT or LONG</td>
<td>The number of rows in each strip (except possibly the last strip, which can contain less rows)</td>
</tr>
<tr class="odd">
<td>Strip Offsets</td>
<td>273</td>
<td>SHORT, LONG or LONG8</td>
<td>For each strip, the byte offset of that strip</td>
</tr>
<tr class="even">
<td>Strip Byte Counts</td>
<td>279</td>
<td>SHORT, LONG or LONG8</td>
<td>For each strip, the number of bytes in that strip</td>
</tr>
</tbody>
</table>

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

<table>
<thead>
<tr class="header">
<th><strong>Field</strong></th>
<th><strong>Tag</strong></th>
<th><strong>Type</strong></th>
<th><strong>Comment</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Image Width</td>
<td>256</td>
<td>SHORT or LONG</td>
<td>The number of columns in the frame</td>
</tr>
<tr class="even">
<td>Image Length</td>
<td>257</td>
<td>SHORT or LONG</td>
<td>The number of rows in the frame</td>
</tr>
<tr class="odd">
<td>Compression</td>
<td>259</td>
<td>SHORT</td>
<td>The compression scheme used; COMPRESSION_NONE (1) for final image data</td>
</tr>
<tr class="even">
<td>Bits per sample</td>
<td>258</td>
<td>SHORT</td>
<td>The number of bits per channel. For final image data, it is 16 bits per sample (2 bytes)</td>
</tr>
<tr class="odd">
<td>Rows Per Strip</td>
<td>278</td>
<td>SHORT or LONG</td>
<td>The number of rows in each strip (except possibly the last strip, which can contain less rows)</td>
</tr>
<tr class="even">
<td>Strip Offsets</td>
<td>273</td>
<td>SHORT, LONG or LONG8</td>
<td>For each strip, the byte offset of that strip</td>
</tr>
<tr class="odd">
<td>Strip Byte Counts</td>
<td>279</td>
<td>SHORT, LONG or LONG8</td>
<td>For each strip, the number of bytes in that strip</td>
</tr>
</tbody>
</table>

Metadata
========

In the first IFD metadata is written that applies to the entire
acquisition. The metadata is written as a tag with ID 65001. This tag
contains an XML string structured as follows:

&lt;metadata&gt;  
&lt;item name=”\[name\]” unit=”\[unit\]”&gt;\[value\]&lt;/item&gt;  
&lt;item name=”\[name2\]”&gt;\[value2\]&lt;/item&gt;  
…  
&lt;/metadata&gt;

The unit attribute may be omitted for values that do not have a unit.

The values present in the field are the following:

<table>
<thead>
<tr class="header">
<th><strong>Name</strong></th>
<th><strong>Unit</strong></th>
<th><strong>Comments</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>numberOfFrames</td>
<td></td>
<td>Total number of frames in the acquisition</td>
</tr>
<tr class="even">
<td>exposureTime</td>
<td>s</td>
<td>Exposure time of the acquisition</td>
</tr>
<tr class="odd">
<td>totalDose</td>
<td>e/pixel</td>
<td>Average amount of electrons per pixel in the entire acquisition</td>
</tr>
<tr class="even">
<td>meanDoseRate</td>
<td>e/pixel/s</td>
<td>Average amount of electrons per pixel per second during the acquisition</td>
</tr>
<tr class="odd">
<td>sensorImageWidth</td>
<td>pixels</td>
<td>Width of the individual images</td>
</tr>
<tr class="even">
<td>sensorImageHeight</td>
<td>pixels</td>
<td>Height of the individual images</td>
</tr>
</tbody>
</table>



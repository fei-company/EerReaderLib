// Copyright (c) 2019 by FEI Company

#ifndef FEIBITSTREAM_H
#define FEIBITSTREAM_H

#include <queue>
#include <iostream>
#include "stdint.h"


typedef uint8_t BitType;
typedef uint16_t MultiBitType;

typedef uint64_t BitStreamWordType;
const int BitStreamWordTypeNBits = 8 * sizeof(BitStreamWordType);


// NB: nWords not used, not checked anyware; responsibility of caller now.
class BitStreamer
{
public:
	BitStreamer(BitStreamWordType* buffer) :
		buffer(buffer), ptr(buffer), bitCounter(BitStreamWordTypeNBits)
	{
		curValue = *ptr;
	}

	// todo this can be optimized.
	MultiBitType getBits(unsigned NBits)
	{
		MultiBitType r = curValue & ((1<<NBits)-1);
		if (bitCounter <= NBits)
		{
			// take a new one.
			curValue = *(++ptr);
			unsigned xtraBitsToRead = NBits - bitCounter; //# Remaining bits.
			r |= ((curValue & ((1<<xtraBitsToRead)-1)) << (bitCounter));
			curValue >>= xtraBitsToRead;
			bitCounter = BitStreamWordTypeNBits - xtraBitsToRead;
		}
		else
		{
			// just shift and done!
			curValue >>= NBits;
			bitCounter -= NBits;
		}
		return r;
	}
 
   void rewind(unsigned NBits)
   {
       //std::cout<<"REWIND "<<bitCounter<<", "<<curValue<<std::endl;
       bitCounter += NBits;
       if (bitCounter > BitStreamWordTypeNBits)
       {
           bitCounter -= BitStreamWordTypeNBits;
           --ptr;
       }
       curValue = (*ptr) >> (BitStreamWordTypeNBits-bitCounter);
       //std::cout<<"/REWINDED "<<bitCounter<<", "<<curValue<<std::endl;
   }

    uint64_t nWordsRead() const
    {
        return (ptr - buffer) + ((bitCounter == BitStreamWordTypeNBits)? 0 : 1);
    }

	BitStreamWordType* buffer;

private:

	BitStreamWordType* ptr;
	BitStreamWordType curValue;
	size_t bitCounter;
};



class BitStreamWriter
{
public:
	BitStreamWriter(BitStreamWordType* buf) :
		ptr(buf), wordPos(0), counter(0)
	{
		*ptr = 0;
    }


	// todo this can be optimized.
	void putBits(MultiBitType v, unsigned NBits)
	{
		v &= ((1 << NBits) - 1); // people are allowed to send words with more bits, for convenience.
        if (wordPos + NBits < BitStreamWordTypeNBits)
        {
            // all fit in current word and 1 bit remains so no need to go to next word pos.
            //put it there, and update word position and done.
            (*ptr) |= (static_cast<BitStreamWordType>(v) << wordPos);
            wordPos += NBits;
        }
        else
        {
            // wordPos+1 most sign.bits still fit in. put them in, update NBits, and go to next word position.
            unsigned NBitsStillFit = BitStreamWordTypeNBits - wordPos;
            (*ptr) |= static_cast<BitStreamWordType>((v & ((1<<NBitsStillFit)-1))) << wordPos;
            v >>= NBitsStillFit;
            // go to next pos;
            ++ptr; ++counter;
            wordPos = NBits - NBitsStillFit;
            *ptr = v;
        }
	}

	size_t getWordSize()
	{
		size_t r = ((size_t)(counter));
		if (wordPos > 0) r+=1;
		return r;
	}

private:
	BitStreamWordType* ptr;
	size_t counter;
    BitStreamWordType wordPos;
};






#endif

/*
 * Copyright (C) 2019 Thermo Fisher Scientific. Do not distribute.
 * This software is provided by Thermo Fisher Scientific to
 * - The Hospital for Sick Children, Toronto, Canada, under confidentiality conditions 
 *   described in the agreed “Equipment Evaluation Agreement” (signed November 2018),
 * - Structura Biotechnology Inc., Toronto, Canada, under confidentiality conditions 
 *   described in the agreed “Mutual Non-Disclosure Agreement” (signed February 2019).
 */

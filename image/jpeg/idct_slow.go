// Copyright 2009 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package jpeg

// This is a Go translation of idct.c from
//
// http://standards.iso.org/ittf/PubliclyAvailableStandards/ISO_IEC_13818-4_2004_Conformance_Testing/Video/verifier/mpeg2decode_960109.tar.gz
//
// which carries the following notice:

/* Copyright (C) 1996, MPEG Software Simulation Group. All Rights Reserved. */

/*
 * Disclaimer of Warranty
 *
 * These software programs are available to the user without any license fee or
 * royalty on an "as is" basis.  The MPEG Software Simulation Group disclaims
 * any and all warranties, whether express, implied, or statuary, including any
 * implied warranties or merchantability or of fitness for a particular
 * purpose.  In no event shall the copyright-holder be liable for any
 * incidental, punitive, or consequential damages of any kind whatsoever
 * arising from the use of these programs.
 *
 * This disclaimer of warranty extends to the user of these programs and user's
 * customers, employees, agents, transferees, successors, and assigns.
 *
 * The MPEG Software Simulation Group does not represent or warrant that the
 * programs furnished hereunder are free of infringement of any third-party
 * patents.
 *
 * Commercial implementations of MPEG-1 and MPEG-2 video, including shareware,
 * are subject to royalty fees to patent holders.  Many of these patents are
 * general enough such that they are unavoidable regardless of implementation
 * design.
 *
 */

// const blockSize = 64 // A DCT block is 8x8.

// type block [blockSize]int32

// idct performs a 2-D Inverse Discrete Cosine Transformation.
//
// The input coefficients should already have been multiplied by the
// appropriate quantization table. We use fixed-point computation, with the
// number of bits for the fractional component varying over the intermediate
// stages.
//
// For more on the actual algorithm, see Z. Wang, "Fast algorithms for the
// discrete W transform and for the discrete Fourier transform", IEEE Trans. on
// ASSP, Vol. ASSP- 32, pp. 803-816, Aug. 1984.

const (
	FIX_0_298631336 = 2446
	FIX_0_390180644 = 3196
	FIX_0_541196100 = 4433
	FIX_0_765366865 = 6270
	FIX_0_899976223 = 7373
	FIX_1_175875602 = 9633
	FIX_1_501321110 = 12299
	FIX_1_847759065 = 15137
	FIX_1_961570560 = 16069
	FIX_2_053119869 = 16819
	FIX_2_562915447 = 20995
	FIX_3_072711026 = 25172
	PASS1_BITS      = 2
	DCTSIZE         = 8
	CONST_BITS      = 13
	ONE             = 1
	RANGE_CENTER    = CenterJSample << 2
	RANGE_MASK      = RANGE_CENTER*2 - 1
)

const (
	MaxJSample     = 255                             // for 8-bit images
	CenterJSample  = 128                             // = 256
	RangeTableSize = RANGE_CENTER*2 + MaxJSample + 1 // full range
)

var range_limit_table = prepareRangeLimitTable()

func range_limit(x int32) int32 {
	return int32(range_limit_table[(CenterJSample+x)&RANGE_MASK])
}

func prepareRangeLimitTable() []uint8 {
	table := make([]uint8, RangeTableSize)

	// Table layout: [0..RangeCenter-1] = 0
	for i := 0; i < RANGE_CENTER; i++ {
		table[i] = 0
	}

	// [RangeCenter..RangeCenter+255] = 0..255
	for i := 0; i <= MaxJSample; i++ {
		table[RANGE_CENTER+i] = uint8(i)
	}

	// [RangeCenter+256..] = 255
	for i := RANGE_CENTER + MaxJSample + 1; i < len(table); i++ {
		table[i] = MaxJSample
	}

	return table
}

func MULTIPLY(a, b int32) int32 {
	return a * b
}

func DEQUANTIZE(a, b int32) int32 {
	return a * b
}

func RIGHT_SHIFT(a, b int32) int32 {
	return a >> b
}

func IRIGHT_SHIFT(a, b int32) int32 {
	return a >> b
}

func idct_slow(src, qt *block) {
	var tmp0, tmp1, tmp2, tmp3 int32
	var tmp10, tmp11, tmp12, tmp13 int32
	var z1, z2, z3 int32
	var workspace [64]int32

	/* Pass 1: process columns from input, store into work array.
	 * Note results are scaled up by sqrt(8) compared to a true IDCT;
	 * furthermore, we scale the results by 2**PASS1_BITS.
	 */

	for ctr := 0; ctr < DCTSIZE; ctr++ {
		inptr := src[ctr:]
		wsptr := workspace[ctr:]
		quantptr := jpegZigzagOrder[ctr:]

		/* Due to quantization, we will usually find that many of the input
		 * coefficients are zero, especially the AC terms.  We can exploit this
		 * by short-circuiting the IDCT calculation for any column in which all
		 * the AC terms are zero.  In that case each output is equal to the
		 * DC coefficient (with scale factor as needed).
		 * With typical images and quantization tables, half or more of the
		 * column DCT calculations can be simplified this way.
		 */
		if inptr[DCTSIZE*1] == 0 && inptr[DCTSIZE*2] == 0 && inptr[DCTSIZE*3] == 0 &&
			inptr[DCTSIZE*4] == 0 && inptr[DCTSIZE*5] == 0 && inptr[DCTSIZE*6] == 0 && inptr[DCTSIZE*7] == 0 {
			/* AC terms all zero */
			dcval := DEQUANTIZE(inptr[DCTSIZE*0], qt[quantptr[DCTSIZE*0]]) << PASS1_BITS

			wsptr[DCTSIZE*0] = dcval
			wsptr[DCTSIZE*1] = dcval
			wsptr[DCTSIZE*2] = dcval
			wsptr[DCTSIZE*3] = dcval
			wsptr[DCTSIZE*4] = dcval
			wsptr[DCTSIZE*5] = dcval
			wsptr[DCTSIZE*6] = dcval
			wsptr[DCTSIZE*7] = dcval

			continue
		}

		/* Even part: reverse the even part of the forward DCT.
		 * The rotator is c(-6).
		 */
		z2 = DEQUANTIZE(inptr[DCTSIZE*0], qt[quantptr[DCTSIZE*0]])
		z3 = DEQUANTIZE(inptr[DCTSIZE*4], qt[quantptr[DCTSIZE*4]])
		z2 <<= CONST_BITS
		z3 <<= CONST_BITS
		/* Add fudge factor here for final descale. */
		z2 += ONE << (CONST_BITS - PASS1_BITS - 1)

		tmp0 = z2 + z3
		tmp1 = z2 - z3

		z2 = DEQUANTIZE(inptr[DCTSIZE*2], qt[quantptr[DCTSIZE*2]])
		z3 = DEQUANTIZE(inptr[DCTSIZE*6], qt[quantptr[DCTSIZE*6]])

		z1 = MULTIPLY(z2+z3, FIX_0_541196100)     /* c6 */
		tmp2 = z1 + MULTIPLY(z2, FIX_0_765366865) /* c2-c6 */
		tmp3 = z1 - MULTIPLY(z3, FIX_1_847759065) /* c2+c6 */

		tmp10 = tmp0 + tmp2
		tmp13 = tmp0 - tmp2
		tmp11 = tmp1 + tmp3
		tmp12 = tmp1 - tmp3

		/* Odd part per figure 8; the matrix is unitary and hence its
		 * transpose is its inverse.  i0..i3 are y7,y5,y3,y1 respectively.
		 */

		tmp0 = DEQUANTIZE(inptr[DCTSIZE*7], qt[quantptr[DCTSIZE*7]])
		tmp1 = DEQUANTIZE(inptr[DCTSIZE*5], qt[quantptr[DCTSIZE*5]])
		tmp2 = DEQUANTIZE(inptr[DCTSIZE*3], qt[quantptr[DCTSIZE*3]])
		tmp3 = DEQUANTIZE(inptr[DCTSIZE*1], qt[quantptr[DCTSIZE*1]])

		z2 = tmp0 + tmp2
		z3 = tmp1 + tmp3

		z1 = MULTIPLY(z2+z3, FIX_1_175875602) /*  c3 */
		z2 = MULTIPLY(z2, -FIX_1_961570560)   /* -c3-c5 */
		z3 = MULTIPLY(z3, -FIX_0_390180644)   /* -c3+c5 */
		z2 += z1
		z3 += z1

		z1 = MULTIPLY(tmp0+tmp3, -FIX_0_899976223) /* -c3+c7 */
		tmp0 = MULTIPLY(tmp0, FIX_0_298631336)     /* -c1+c3+c5-c7 */
		tmp3 = MULTIPLY(tmp3, FIX_1_501321110)     /*  c1+c3-c5-c7 */
		tmp0 += z1 + z2
		tmp3 += z1 + z3

		z1 = MULTIPLY(tmp1+tmp2, -FIX_2_562915447) /* -c1-c3 */
		tmp1 = MULTIPLY(tmp1, FIX_2_053119869)     /*  c1+c3-c5+c7 */
		tmp2 = MULTIPLY(tmp2, FIX_3_072711026)     /*  c1+c3+c5-c7 */
		tmp1 += z1 + z3
		tmp2 += z1 + z2

		/* Final output stage: inputs are tmp10..tmp13, tmp0..tmp3 */

		wsptr[DCTSIZE*0] = RIGHT_SHIFT(tmp10+tmp3, CONST_BITS-PASS1_BITS)
		wsptr[DCTSIZE*7] = RIGHT_SHIFT(tmp10-tmp3, CONST_BITS-PASS1_BITS)
		wsptr[DCTSIZE*1] = RIGHT_SHIFT(tmp11+tmp2, CONST_BITS-PASS1_BITS)
		wsptr[DCTSIZE*6] = RIGHT_SHIFT(tmp11-tmp2, CONST_BITS-PASS1_BITS)
		wsptr[DCTSIZE*2] = RIGHT_SHIFT(tmp12+tmp1, CONST_BITS-PASS1_BITS)
		wsptr[DCTSIZE*5] = RIGHT_SHIFT(tmp12-tmp1, CONST_BITS-PASS1_BITS)
		wsptr[DCTSIZE*3] = RIGHT_SHIFT(tmp13+tmp0, CONST_BITS-PASS1_BITS)
		wsptr[DCTSIZE*4] = RIGHT_SHIFT(tmp13-tmp0, CONST_BITS-PASS1_BITS)
	}

	/* Pass 2: process rows from work array, store into output array.
	 * Note that we must descale the results by a factor of 8 == 2**3,
	 * and also undo the PASS1_BITS scaling.
	 */
	for ctr := 0; ctr < DCTSIZE; ctr++ {
		wsptr := workspace[ctr*DCTSIZE:]
		outptr := src[ctr*DCTSIZE:]

		/* Add range center and fudge factor for final descale and range-limit. */
		z2 = wsptr[0] +
			(((RANGE_CENTER) << (PASS1_BITS + 3)) +
				(ONE << (PASS1_BITS + 2)))

		/* Rows of zeroes can be exploited in the same way as we did with columns.
		 * However, the column calculation has created many nonzero AC terms, so
		 * the simplification applies less often (typically 5% to 10% of the time).
		 * On machines with very fast multiplication, it's possible that the
		 * test takes more time than it's worth.  In that case this section
		 * may be commented out.
		 */

		if wsptr[1] == 0 && wsptr[2] == 0 && wsptr[3] == 0 && wsptr[4] == 0 &&
			wsptr[5] == 0 && wsptr[6] == 0 && wsptr[7] == 0 {
			dcval := range_limit(RIGHT_SHIFT(z2, PASS1_BITS+3))
			outptr[0] = dcval
			outptr[1] = dcval
			outptr[2] = dcval
			outptr[3] = dcval
			outptr[4] = dcval
			outptr[5] = dcval
			outptr[6] = dcval
			outptr[7] = dcval
			continue
		}
		/* Even part: reverse the even part of the forward DCT.
		 * The rotator is c(-6).
		 */

		z3 = wsptr[4]

		tmp0 = (z2 + z3) << CONST_BITS
		tmp1 = (z2 - z3) << CONST_BITS

		z2 = wsptr[2]
		z3 = wsptr[6]

		z1 = MULTIPLY(z2+z3, FIX_0_541196100)     /* c6 */
		tmp2 = z1 + MULTIPLY(z2, FIX_0_765366865) /* c2-c6 */
		tmp3 = z1 - MULTIPLY(z3, FIX_1_847759065) /* c2+c6 */

		tmp10 = tmp0 + tmp2
		tmp13 = tmp0 - tmp2
		tmp11 = tmp1 + tmp3
		tmp12 = tmp1 - tmp3

		/* Odd part per figure 8; the matrix is unitary and hence its
		 * transpose is its inverse.  i0..i3 are y7,y5,y3,y1 respectively.
		 */

		tmp0 = wsptr[7]
		tmp1 = wsptr[5]
		tmp2 = wsptr[3]
		tmp3 = wsptr[1]

		z2 = tmp0 + tmp2
		z3 = tmp1 + tmp3

		z1 = MULTIPLY(z2+z3, FIX_1_175875602) /*  c3 */
		z2 = MULTIPLY(z2, -FIX_1_961570560)   /* -c3-c5 */
		z3 = MULTIPLY(z3, -FIX_0_390180644)   /* -c3+c5 */
		z2 += z1
		z3 += z1

		z1 = MULTIPLY(tmp0+tmp3, -FIX_0_899976223) /* -c3+c7 */
		tmp0 = MULTIPLY(tmp0, FIX_0_298631336)     /* -c1+c3+c5-c7 */
		tmp3 = MULTIPLY(tmp3, FIX_1_501321110)     /*  c1+c3-c5-c7 */
		tmp0 += z1 + z2
		tmp3 += z1 + z3

		z1 = MULTIPLY(tmp1+tmp2, -FIX_2_562915447) /* -c1-c3 */
		tmp1 = MULTIPLY(tmp1, FIX_2_053119869)     /*  c1+c3-c5+c7 */
		tmp2 = MULTIPLY(tmp2, FIX_3_072711026)     /*  c1+c3+c5-c7 */
		tmp1 += z1 + z3
		tmp2 += z1 + z2

		/* Final output stage: inputs are tmp10..tmp13, tmp0..tmp3 */

		outptr[0] = range_limit(RIGHT_SHIFT(tmp10+tmp3, CONST_BITS+PASS1_BITS+3))
		outptr[7] = range_limit(RIGHT_SHIFT(tmp10-tmp3, CONST_BITS+PASS1_BITS+3))
		outptr[1] = range_limit(RIGHT_SHIFT(tmp11+tmp2, CONST_BITS+PASS1_BITS+3))
		outptr[6] = range_limit(RIGHT_SHIFT(tmp11-tmp2, CONST_BITS+PASS1_BITS+3))
		outptr[2] = range_limit(RIGHT_SHIFT(tmp12+tmp1, CONST_BITS+PASS1_BITS+3))
		outptr[5] = range_limit(RIGHT_SHIFT(tmp12-tmp1, CONST_BITS+PASS1_BITS+3))
		outptr[3] = range_limit(RIGHT_SHIFT(tmp13+tmp0, CONST_BITS+PASS1_BITS+3))
		outptr[4] = range_limit(RIGHT_SHIFT(tmp13-tmp0, CONST_BITS+PASS1_BITS+3))
	}
}

var jpegZigzagOrder = [blockSize]int{
	0, 1, 5, 6, 14, 15, 27, 28,
	2, 4, 7, 13, 16, 26, 29, 42,
	3, 8, 12, 17, 25, 30, 41, 43,
	9, 11, 18, 24, 31, 40, 44, 53,
	10, 19, 23, 32, 39, 45, 52, 54,
	20, 22, 33, 38, 46, 51, 55, 60,
	21, 34, 37, 47, 50, 56, 59, 61,
	35, 36, 48, 49, 57, 58, 62, 63,
}

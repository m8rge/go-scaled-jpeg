package jpeg

// This is a Go translation of jidctint.c from
//
// https://www.ijg.org/files/jpegsrc.v9f.tar.gz
//
// which carries the following notice:

/*
 * jidctint.c
 *
 * Copyright (C) 1991-1998, Thomas G. Lane.
 * Modification developed 2002-2018 by Guido Vollbeding.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains a slow-but-accurate integer implementation of the
 * inverse DCT (Discrete Cosine Transform).  In the IJG code, this routine
 * must also perform dequantization of the input coefficients.
 *
 * A 2-D IDCT can be done by 1-D IDCT on each column followed by 1-D IDCT
 * on each row (or vice versa, but it's more convenient to emit a row at
 * a time).  Direct algorithms are also available, but they are much more
 * complex and seem not to be any faster when reduced to code.
 *
 * This implementation is based on an algorithm described in
 *   C. Loeffler, A. Ligtenberg and G. Moschytz, "Practical Fast 1-D DCT
 *   Algorithms with 11 Multiplications", Proc. Int'l. Conf. on Acoustics,
 *   Speech, and Signal Processing 1989 (ICASSP '89), pp. 988-991.
 * The primary algorithm described there uses 11 multiplies and 29 adds.
 * We use their alternate method with 12 multiplies and 32 adds.
 * The advantage of this method is that no data path contains more than one
 * multiplication; this allows a very simple and accurate implementation in
 * scaled fixed-point arithmetic, with a minimal number of shifts.
 *
 * We also provide IDCT routines with various output sample block sizes for
 * direct resolution reduction or enlargement and for direct resolving the
 * common 2x1 and 1x2 subsampling cases without additional resampling: NxN
 * (N=1...16), 2NxN, and Nx2N (N=1...8) pixels for one 8x8 input DCT block.
 *
 * For N<8 we simply take the corresponding low-frequency coefficients of
 * the 8x8 input DCT block and apply an NxN point IDCT on the sub-block
 * to yield the downscaled outputs.
 * This can be seen as direct low-pass downsampling from the DCT domain
 * point of view rather than the usual spatial domain point of view,
 * yielding significant computational savings and results at least
 * as good as common bilinear (averaging) spatial downsampling.
 *
 * For N>8 we apply a partial NxN IDCT on the 8 input coefficients as
 * lower frequencies and higher frequencies assumed to be zero.
 * It turns out that the computational effort is similar to the 8x8 IDCT
 * regarding the output size.
 * Furthermore, the scaling and descaling is the same for all IDCT sizes.
 *
 * CAUTION: We rely on the FIX() macro except for the N=1,2,4,8 cases
 * since there would be too many additional constants to pre-calculate.
 */

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
	RANGE_CENTER    = CENTERJSAMPLE << 2
	RANGE_MASK      = RANGE_CENTER*2 - 1
	MAXJSAMPLE      = 255 // for 8-bit images
	CENTERJSAMPLE   = 128
	RangeTableSize  = RANGE_CENTER*2 + MAXJSAMPLE + 1 // full range
)

var range_limit_table = prepareRangeLimitTable()

func range_limit(x int32) int32 {
	return int32(range_limit_table[(CENTERJSAMPLE+x)&RANGE_MASK])
}

func prepareRangeLimitTable() []uint8 {
	table := make([]uint8, RangeTableSize)

	// Table layout: [0..RangeCenter-1] = 0
	for i := 0; i < RANGE_CENTER; i++ {
		table[i] = 0
	}

	// [RangeCenter..RangeCenter+255] = 0..255
	for i := 0; i <= MAXJSAMPLE; i++ {
		table[RANGE_CENTER+i] = uint8(i)
	}

	// [RangeCenter+256..] = 255
	for i := RANGE_CENTER + MAXJSAMPLE + 1; i < len(table); i++ {
		table[i] = MAXJSAMPLE
	}

	return table
}

func MULTIPLY(a, b int32) int32 {
	return a * b
}

/* Dequantize a coefficient by multiplying it by the multiplier-table
 * entry; produce an int result.  In this module, both inputs and result
 * are 16 bits or less, so either int or short multiply will work.
 */
func DEQUANTIZE(a, b int32) int32 {
	return a * b
}

func RIGHT_SHIFT(a, b int32) int32 {
	return a >> b
}

func IRIGHT_SHIFT(a, b int32) int32 {
	return a >> b
}

const CONST_SCALE = ONE << CONST_BITS

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

/* Convert a positive real constant to an integer scaled by CONST_SCALE.
 * Caution: some C compilers fail to reduce "FIX(constant)" at compile time,
 * thus causing a lot of useless floating-point operations at run time.
 */
func FIX(a float32) int32 {
	return int32(a*CONST_SCALE + 0.5)
}

/*
 * Perform dequantization and inverse DCT on one block of coefficients.
 *
 * Optimized algorithm with 12 multiplications in the 1-D kernel.
 * cK represents sqrt(2) * cos(K*pi/16).
 */
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

/*
 * Perform dequantization and inverse DCT on one block of coefficients,
 * producing a reduced-size 7x7 output block.
 *
 * Optimized algorithm with 12 multiplications in the 1-D kernel.
 * cK represents sqrt(2) * cos(K*pi/14).
 */
func jpeg_idct_7x7(src, qt *block) {
	var tmp0, tmp1, tmp2, tmp10, tmp11, tmp12, tmp13 int32
	var z1, z2, z3 int32
	var workspace [7 * 7]int32 /* buffers data between passes * int32

	/* Pass 1: process columns from input, store into work array. */

	for ctr := 0; ctr < 7; ctr++ {
		inptr := src[ctr:]
		wsptr := workspace[ctr:]
		quantptr := jpegZigzagOrder[ctr:]

		/* Even part */

		tmp13 = DEQUANTIZE(inptr[DCTSIZE*0], qt[quantptr[DCTSIZE*0]])
		tmp13 <<= CONST_BITS
		/* Add fudge factor here for final descale. */
		tmp13 += ONE << (CONST_BITS - PASS1_BITS - 1)

		z1 = DEQUANTIZE(inptr[DCTSIZE*2], qt[quantptr[DCTSIZE*2]])
		z2 = DEQUANTIZE(inptr[DCTSIZE*4], qt[quantptr[DCTSIZE*4]])
		z3 = DEQUANTIZE(inptr[DCTSIZE*6], qt[quantptr[DCTSIZE*6]])

		tmp10 = MULTIPLY(z2-z3, FIX(0.881747734))                      /* c4 */
		tmp12 = MULTIPLY(z1-z2, FIX(0.314692123))                      /* c6 */
		tmp11 = tmp10 + tmp12 + tmp13 - MULTIPLY(z2, FIX(1.841218003)) /* c2+c4-c6 */
		tmp0 = z1 + z3
		z2 -= tmp0
		tmp0 = MULTIPLY(tmp0, FIX(1.274162392)) + tmp13 /* c2 */
		tmp10 += tmp0 - MULTIPLY(z3, FIX(0.077722536))  /* c2-c4-c6 */
		tmp12 += tmp0 - MULTIPLY(z1, FIX(2.470602249))  /* c2+c4+c6 */
		tmp13 += MULTIPLY(z2, FIX(1.414213562))         /* c0 */

		/* Odd part */

		z1 = DEQUANTIZE(inptr[DCTSIZE*1], qt[quantptr[DCTSIZE*1]])
		z2 = DEQUANTIZE(inptr[DCTSIZE*3], qt[quantptr[DCTSIZE*3]])
		z3 = DEQUANTIZE(inptr[DCTSIZE*5], qt[quantptr[DCTSIZE*5]])

		tmp1 = MULTIPLY(z1+z2, FIX(0.935414347)) /* (c3+c1-c5)/2 */
		tmp2 = MULTIPLY(z1-z2, FIX(0.170262339)) /* (c3+c5-c1)/2 */
		tmp0 = tmp1 - tmp2
		tmp1 += tmp2
		tmp2 = MULTIPLY(z2+z3, -FIX(1.378756276)) /* -c1 */
		tmp1 += tmp2
		z2 = MULTIPLY(z1+z3, FIX(0.613604268)) /* c5 */
		tmp0 += z2
		tmp2 += z2 + MULTIPLY(z3, FIX(1.870828693)) /* c3+c1-c5 */

		/* Final output stage */

		wsptr[7*0] = RIGHT_SHIFT(tmp10+tmp0, CONST_BITS-PASS1_BITS)
		wsptr[7*6] = RIGHT_SHIFT(tmp10-tmp0, CONST_BITS-PASS1_BITS)
		wsptr[7*1] = RIGHT_SHIFT(tmp11+tmp1, CONST_BITS-PASS1_BITS)
		wsptr[7*5] = RIGHT_SHIFT(tmp11-tmp1, CONST_BITS-PASS1_BITS)
		wsptr[7*2] = RIGHT_SHIFT(tmp12+tmp2, CONST_BITS-PASS1_BITS)
		wsptr[7*4] = RIGHT_SHIFT(tmp12-tmp2, CONST_BITS-PASS1_BITS)
		wsptr[7*3] = RIGHT_SHIFT(tmp13, CONST_BITS-PASS1_BITS)
	}

	/* Pass 2: process 7 rows from work array, store into output array. */

	for ctr := 0; ctr < 7; ctr++ {
		wsptr := workspace[ctr*7:]
		outptr := src[ctr*7:]

		/* Even part */

		/* Add range center and fudge factor for final descale and range-limit. */
		tmp13 = wsptr[0] +
			(((RANGE_CENTER) << (PASS1_BITS + 3)) +
				(ONE << (PASS1_BITS + 2)))
		tmp13 <<= CONST_BITS

		z1 = wsptr[2]
		z2 = wsptr[4]
		z3 = wsptr[6]

		tmp10 = MULTIPLY(z2-z3, FIX(0.881747734))                      /* c4 */
		tmp12 = MULTIPLY(z1-z2, FIX(0.314692123))                      /* c6 */
		tmp11 = tmp10 + tmp12 + tmp13 - MULTIPLY(z2, FIX(1.841218003)) /* c2+c4-c6 */
		tmp0 = z1 + z3
		z2 -= tmp0
		tmp0 = MULTIPLY(tmp0, FIX(1.274162392)) + tmp13 /* c2 */
		tmp10 += tmp0 - MULTIPLY(z3, FIX(0.077722536))  /* c2-c4-c6 */
		tmp12 += tmp0 - MULTIPLY(z1, FIX(2.470602249))  /* c2+c4+c6 */
		tmp13 += MULTIPLY(z2, FIX(1.414213562))         /* c0 */

		/* Odd part */

		z1 = wsptr[1]
		z2 = wsptr[3]
		z3 = wsptr[5]

		tmp1 = MULTIPLY(z1+z2, FIX(0.935414347)) /* (c3+c1-c5)/2 */
		tmp2 = MULTIPLY(z1-z2, FIX(0.170262339)) /* (c3+c5-c1)/2 */
		tmp0 = tmp1 - tmp2
		tmp1 += tmp2
		tmp2 = MULTIPLY(z2+z3, -FIX(1.378756276)) /* -c1 */
		tmp1 += tmp2
		z2 = MULTIPLY(z1+z3, FIX(0.613604268)) /* c5 */
		tmp0 += z2
		tmp2 += z2 + MULTIPLY(z3, FIX(1.870828693)) /* c3+c1-c5 */

		/* Final output stage */

		outptr[0] = range_limit(RIGHT_SHIFT(tmp10+tmp0, CONST_BITS+PASS1_BITS+3))
		outptr[6] = range_limit(RIGHT_SHIFT(tmp10-tmp0, CONST_BITS+PASS1_BITS+3))
		outptr[1] = range_limit(RIGHT_SHIFT(tmp11+tmp1, CONST_BITS+PASS1_BITS+3))
		outptr[5] = range_limit(RIGHT_SHIFT(tmp11-tmp1, CONST_BITS+PASS1_BITS+3))
		outptr[2] = range_limit(RIGHT_SHIFT(tmp12+tmp2, CONST_BITS+PASS1_BITS+3))
		outptr[4] = range_limit(RIGHT_SHIFT(tmp12-tmp2, CONST_BITS+PASS1_BITS+3))
		outptr[3] = range_limit(RIGHT_SHIFT(tmp13, CONST_BITS+PASS1_BITS+3))

	}
}

/*
 * Perform dequantization and inverse DCT on one block of coefficients,
 * producing a reduced-size 6x6 output block.
 *
 * Optimized algorithm with 3 multiplications in the 1-D kernel.
 * cK represents sqrt(2) * cos(K*pi/12).
 */

func jpeg_idct_6x6(src, qt *block) {
	var tmp0, tmp1, tmp2, tmp10, tmp11, tmp12 int32
	var z1, z2, z3 int32

	var workspace [6 * 6]int32 /* buffers data between passes */

	/* Pass 1: process columns from input, store into work array. */

	for ctr := 0; ctr < 6; ctr++ {
		inptr := src[ctr:]
		wsptr := workspace[ctr:]
		quantptr := jpegZigzagOrder[ctr:]
		/* Even part */

		tmp0 = DEQUANTIZE(inptr[DCTSIZE*0], qt[quantptr[DCTSIZE*0]])
		tmp0 <<= CONST_BITS
		/* Add fudge factor here for final descale. */
		tmp0 += ONE << (CONST_BITS - PASS1_BITS - 1)
		tmp2 = DEQUANTIZE(inptr[DCTSIZE*4], qt[quantptr[DCTSIZE*4]])
		tmp10 = MULTIPLY(tmp2, FIX(0.707106781)) /* c4 */
		tmp1 = tmp0 + tmp10
		tmp11 = RIGHT_SHIFT(tmp0-tmp10-tmp10, CONST_BITS-PASS1_BITS)
		tmp10 = DEQUANTIZE(inptr[DCTSIZE*2], qt[quantptr[DCTSIZE*2]])
		tmp0 = MULTIPLY(tmp10, FIX(1.224744871)) /* c2 */
		tmp10 = tmp1 + tmp0
		tmp12 = tmp1 - tmp0

		/* Odd part */

		z1 = DEQUANTIZE(inptr[DCTSIZE*1], qt[quantptr[DCTSIZE*1]])
		z2 = DEQUANTIZE(inptr[DCTSIZE*3], qt[quantptr[DCTSIZE*3]])
		z3 = DEQUANTIZE(inptr[DCTSIZE*5], qt[quantptr[DCTSIZE*5]])
		tmp1 = MULTIPLY(z1+z3, FIX(0.366025404)) /* c5 */
		tmp0 = tmp1 + ((z1 + z2) << CONST_BITS)
		tmp2 = tmp1 + ((z3 - z2) << CONST_BITS)
		tmp1 = (z1 - z2 - z3) << PASS1_BITS

		/* Final output stage */

		wsptr[6*0] = RIGHT_SHIFT(tmp10+tmp0, CONST_BITS-PASS1_BITS)
		wsptr[6*5] = RIGHT_SHIFT(tmp10-tmp0, CONST_BITS-PASS1_BITS)
		wsptr[6*1] = tmp11 + tmp1
		wsptr[6*4] = tmp11 - tmp1
		wsptr[6*2] = RIGHT_SHIFT(tmp12+tmp2, CONST_BITS-PASS1_BITS)
		wsptr[6*3] = RIGHT_SHIFT(tmp12-tmp2, CONST_BITS-PASS1_BITS)
	}

	/* Pass 2: process 6 rows from work array, store into output array. */

	for ctr := 0; ctr < 6; ctr++ {
		wsptr := workspace[ctr*6:]
		outptr := src[ctr*6:]

		/* Even part */

		/* Add range center and fudge factor for final descale and range-limit. */
		tmp0 = wsptr[0] +
			(((RANGE_CENTER) << (PASS1_BITS + 3)) +
				(ONE << (PASS1_BITS + 2)))
		tmp0 <<= CONST_BITS
		tmp2 = wsptr[4]
		tmp10 = MULTIPLY(tmp2, FIX(0.707106781)) /* c4 */
		tmp1 = tmp0 + tmp10
		tmp11 = tmp0 - tmp10 - tmp10
		tmp10 = wsptr[2]
		tmp0 = MULTIPLY(tmp10, FIX(1.224744871)) /* c2 */
		tmp10 = tmp1 + tmp0
		tmp12 = tmp1 - tmp0

		/* Odd part */

		z1 = wsptr[1]
		z2 = wsptr[3]
		z3 = wsptr[5]
		tmp1 = MULTIPLY(z1+z3, FIX(0.366025404)) /* c5 */
		tmp0 = tmp1 + ((z1 + z2) << CONST_BITS)
		tmp2 = tmp1 + ((z3 - z2) << CONST_BITS)
		tmp1 = (z1 - z2 - z3) << CONST_BITS

		/* Final output stage */

		outptr[0] = range_limit(RIGHT_SHIFT(tmp10+tmp0, CONST_BITS+PASS1_BITS+3))
		outptr[5] = range_limit(RIGHT_SHIFT(tmp10-tmp0, CONST_BITS+PASS1_BITS+3))
		outptr[1] = range_limit(RIGHT_SHIFT(tmp11+tmp1, CONST_BITS+PASS1_BITS+3))
		outptr[4] = range_limit(RIGHT_SHIFT(tmp11-tmp1, CONST_BITS+PASS1_BITS+3))
		outptr[2] = range_limit(RIGHT_SHIFT(tmp12+tmp2, CONST_BITS+PASS1_BITS+3))
		outptr[3] = range_limit(RIGHT_SHIFT(tmp12-tmp2, CONST_BITS+PASS1_BITS+3))

	}
}

/*
 * Perform dequantization and inverse DCT on one block of coefficients,
 * producing a reduced-size 5x5 output block.
 *
 * Optimized algorithm with 5 multiplications in the 1-D kernel.
 * cK represents sqrt(2) * cos(K*pi/10).
 */

func jpeg_idct_5x5(src, qt *block) {
	var tmp0, tmp1, tmp10, tmp11, tmp12 int32
	var z1, z2, z3 int32

	var workspace [5 * 5]int32 /* buffers data between passes */

	/* Pass 1: process columns from input, store into work array. */

	for ctr := 0; ctr < 5; ctr++ {
		inptr := src[ctr:]
		wsptr := workspace[ctr:]
		quantptr := jpegZigzagOrder[ctr:]
		/* Even part */

		tmp12 = DEQUANTIZE(inptr[DCTSIZE*0], qt[quantptr[DCTSIZE*0]])
		tmp12 <<= CONST_BITS
		/* Add fudge factor here for final descale. */
		tmp12 += ONE << (CONST_BITS - PASS1_BITS - 1)
		tmp0 = DEQUANTIZE(inptr[DCTSIZE*2], qt[quantptr[DCTSIZE*2]])
		tmp1 = DEQUANTIZE(inptr[DCTSIZE*4], qt[quantptr[DCTSIZE*4]])
		z1 = MULTIPLY(tmp0+tmp1, FIX(0.790569415)) /* (c2+c4)/2 */
		z2 = MULTIPLY(tmp0-tmp1, FIX(0.353553391)) /* (c2-c4)/2 */
		z3 = tmp12 + z2
		tmp10 = z3 + z1
		tmp11 = z3 - z1
		tmp12 -= z2 << 2

		/* Odd part */

		z2 = DEQUANTIZE(inptr[DCTSIZE*1], qt[quantptr[DCTSIZE*1]])
		z3 = DEQUANTIZE(inptr[DCTSIZE*3], qt[quantptr[DCTSIZE*3]])

		z1 = MULTIPLY(z2+z3, FIX(0.831253876))     /* c3 */
		tmp0 = z1 + MULTIPLY(z2, FIX(0.513743148)) /* c1-c3 */
		tmp1 = z1 - MULTIPLY(z3, FIX(2.176250899)) /* c1+c3 */

		/* Final output stage */

		wsptr[5*0] = RIGHT_SHIFT(tmp10+tmp0, CONST_BITS-PASS1_BITS)
		wsptr[5*4] = RIGHT_SHIFT(tmp10-tmp0, CONST_BITS-PASS1_BITS)
		wsptr[5*1] = RIGHT_SHIFT(tmp11+tmp1, CONST_BITS-PASS1_BITS)
		wsptr[5*3] = RIGHT_SHIFT(tmp11-tmp1, CONST_BITS-PASS1_BITS)
		wsptr[5*2] = RIGHT_SHIFT(tmp12, CONST_BITS-PASS1_BITS)
	}

	/* Pass 2: process 5 rows from work array, store into output array. */

	for ctr := 0; ctr < 5; ctr++ {
		wsptr := workspace[ctr*5:]
		outptr := src[ctr*5:]

		/* Even part */

		/* Add range center and fudge factor for final descale and range-limit. */
		tmp12 = wsptr[0] +
			(((RANGE_CENTER) << (PASS1_BITS + 3)) +
				(ONE << (PASS1_BITS + 2)))
		tmp12 <<= CONST_BITS
		tmp0 = wsptr[2]
		tmp1 = wsptr[4]
		z1 = MULTIPLY(tmp0+tmp1, FIX(0.790569415)) /* (c2+c4)/2 */
		z2 = MULTIPLY(tmp0-tmp1, FIX(0.353553391)) /* (c2-c4)/2 */
		z3 = tmp12 + z2
		tmp10 = z3 + z1
		tmp11 = z3 - z1
		tmp12 -= z2 << 2

		/* Odd part */

		z2 = wsptr[1]
		z3 = wsptr[3]

		z1 = MULTIPLY(z2+z3, FIX(0.831253876))     /* c3 */
		tmp0 = z1 + MULTIPLY(z2, FIX(0.513743148)) /* c1-c3 */
		tmp1 = z1 - MULTIPLY(z3, FIX(2.176250899)) /* c1+c3 */

		/* Final output stage */

		outptr[0] = range_limit(RIGHT_SHIFT(tmp10+tmp0, CONST_BITS+PASS1_BITS+3))
		outptr[4] = range_limit(RIGHT_SHIFT(tmp10-tmp0, CONST_BITS+PASS1_BITS+3))
		outptr[1] = range_limit(RIGHT_SHIFT(tmp11+tmp1, CONST_BITS+PASS1_BITS+3))
		outptr[3] = range_limit(RIGHT_SHIFT(tmp11-tmp1, CONST_BITS+PASS1_BITS+3))
		outptr[2] = range_limit(RIGHT_SHIFT(tmp12, CONST_BITS+PASS1_BITS+3))

	}
}

func jpeg_idct_4x4(src, qt *block) {
	var tmp0, tmp2, tmp10, tmp12 int32
	var z1, z2, z3 int32

	var workspace [4 * 4]int32 /* buffers data between passes */

	/* Pass 1: process columns from input, store into work array. */

	for ctr := 0; ctr < 4; ctr++ {
		inptr := src[ctr:]
		wsptr := workspace[ctr:]
		quantptr := jpegZigzagOrder[ctr:]
		/* Even part */

		tmp0 = DEQUANTIZE(inptr[DCTSIZE*0], qt[quantptr[DCTSIZE*0]])
		tmp2 = DEQUANTIZE(inptr[DCTSIZE*2], qt[quantptr[DCTSIZE*2]])

		tmp10 = (tmp0 + tmp2) << PASS1_BITS
		tmp12 = (tmp0 - tmp2) << PASS1_BITS

		/* Odd part */
		/* Same rotation as in the even part of the 8x8 LL&M IDCT */

		z2 = DEQUANTIZE(inptr[DCTSIZE*1], qt[quantptr[DCTSIZE*1]])
		z3 = DEQUANTIZE(inptr[DCTSIZE*3], qt[quantptr[DCTSIZE*3]])

		z1 = MULTIPLY(z2+z3, FIX_0_541196100) /* c6 */
		/* Add fudge factor here for final descale. */
		z1 += ONE << (CONST_BITS - PASS1_BITS - 1)
		tmp0 = RIGHT_SHIFT(z1+MULTIPLY(z2, FIX_0_765366865), /* c2-c6 */
			CONST_BITS-PASS1_BITS)
		tmp2 = RIGHT_SHIFT(z1-MULTIPLY(z3, FIX_1_847759065), /* c2+c6 */
			CONST_BITS-PASS1_BITS)

		/* Final output stage */

		wsptr[4*0] = tmp10 + tmp0
		wsptr[4*3] = tmp10 - tmp0
		wsptr[4*1] = tmp12 + tmp2
		wsptr[4*2] = tmp12 - tmp2
	}

	/* Pass 2: process 4 rows from work array, store into output array. */

	for ctr := 0; ctr < 4; ctr++ {
		wsptr := workspace[ctr*4:]
		outptr := src[ctr*4:]

		/* Even part */

		/* Add range center and fudge factor for final descale and range-limit. */
		tmp0 = wsptr[0] +
			(((RANGE_CENTER) << (PASS1_BITS + 3)) +
				(ONE << (PASS1_BITS + 2)))
		tmp2 = wsptr[2]

		tmp10 = (tmp0 + tmp2) << CONST_BITS
		tmp12 = (tmp0 - tmp2) << CONST_BITS

		/* Odd part */
		/* Same rotation as in the even part of the 8x8 LL&M IDCT */

		z2 = wsptr[1]
		z3 = wsptr[3]

		z1 = MULTIPLY(z2+z3, FIX_0_541196100)     /* c6 */
		tmp0 = z1 + MULTIPLY(z2, FIX_0_765366865) /* c2-c6 */
		tmp2 = z1 - MULTIPLY(z3, FIX_1_847759065) /* c2+c6 */

		/* Final output stage */

		outptr[0] = range_limit(RIGHT_SHIFT(tmp10+tmp0, CONST_BITS+PASS1_BITS+3))
		outptr[3] = range_limit(RIGHT_SHIFT(tmp10-tmp0, CONST_BITS+PASS1_BITS+3))
		outptr[1] = range_limit(RIGHT_SHIFT(tmp12+tmp2, CONST_BITS+PASS1_BITS+3))
		outptr[2] = range_limit(RIGHT_SHIFT(tmp12-tmp2, CONST_BITS+PASS1_BITS+3))

	}
}

/*
 * Perform dequantization and inverse DCT on one block of coefficients,
 * producing a reduced-size 3x3 output block.
 *
 * Optimized algorithm with 2 multiplications in the 1-D kernel.
 * cK represents sqrt(2) * cos(K*pi/6).
 */

func jpeg_idct_3x3(src, qt *block) {
	var tmp0, tmp2, tmp10, tmp12 int32

	var workspace [3 * 3]int32 /* buffers data between passes */

	/* Pass 1: process columns from input, store into work array. */

	for ctr := 0; ctr < 3; ctr++ {
		inptr := src[ctr:]
		wsptr := workspace[ctr:]
		quantptr := jpegZigzagOrder[ctr:]
		/* Even part */

		tmp0 = DEQUANTIZE(inptr[DCTSIZE*0], qt[quantptr[DCTSIZE*0]])
		tmp0 <<= CONST_BITS
		/* Add fudge factor here for final descale. */
		tmp0 += ONE << (CONST_BITS - PASS1_BITS - 1)
		tmp2 = DEQUANTIZE(inptr[DCTSIZE*2], qt[quantptr[DCTSIZE*2]])
		tmp12 = MULTIPLY(tmp2, FIX(0.707106781)) /* c2 */
		tmp10 = tmp0 + tmp12
		tmp2 = tmp0 - tmp12 - tmp12

		/* Odd part */

		tmp12 = DEQUANTIZE(inptr[DCTSIZE*1], qt[quantptr[DCTSIZE*1]])
		tmp0 = MULTIPLY(tmp12, FIX(1.224744871)) /* c1 */

		/* Final output stage */

		wsptr[3*0] = RIGHT_SHIFT(tmp10+tmp0, CONST_BITS-PASS1_BITS)
		wsptr[3*2] = RIGHT_SHIFT(tmp10-tmp0, CONST_BITS-PASS1_BITS)
		wsptr[3*1] = RIGHT_SHIFT(tmp2, CONST_BITS-PASS1_BITS)
	}

	/* Pass 2: process 3 rows from work array, store into output array. */

	for ctr := 0; ctr < 3; ctr++ {
		wsptr := workspace[ctr*3:]
		outptr := src[ctr*3:]

		/* Even part */

		/* Add range center and fudge factor for final descale and range-limit. */
		tmp0 = wsptr[0] +
			(((RANGE_CENTER) << (PASS1_BITS + 3)) +
				(ONE << (PASS1_BITS + 2)))
		tmp0 <<= CONST_BITS
		tmp2 = wsptr[2]
		tmp12 = MULTIPLY(tmp2, FIX(0.707106781)) /* c2 */
		tmp10 = tmp0 + tmp12
		tmp2 = tmp0 - tmp12 - tmp12

		/* Odd part */

		tmp12 = wsptr[1]
		tmp0 = MULTIPLY(tmp12, FIX(1.224744871)) /* c1 */

		/* Final output stage */

		outptr[0] = range_limit(RIGHT_SHIFT(tmp10+tmp0, CONST_BITS+PASS1_BITS+3))
		outptr[2] = range_limit(RIGHT_SHIFT(tmp10-tmp0, CONST_BITS+PASS1_BITS+3))
		outptr[1] = range_limit(RIGHT_SHIFT(tmp2, CONST_BITS+PASS1_BITS+3))

	}
}

/*
 * Perform dequantization and inverse DCT on one block of coefficients,
 * producing a reduced-size 2x2 output block.
 *
 * Multiplication-less algorithm.
 */

func jpeg_idct_2x2(src, qt *block) {
	var tmp0, tmp1, tmp2, tmp3, tmp4, tmp5 int32

	coef_block := src[:]
	quantptr := jpegZigzagOrder[:]

	/* Pass 1: process columns from input. */

	/* Column 0 */
	tmp4 = DEQUANTIZE(coef_block[DCTSIZE*0], qt[quantptr[DCTSIZE*0]])
	tmp5 = DEQUANTIZE(coef_block[DCTSIZE*1], qt[quantptr[DCTSIZE*1]])
	/* Add range center and fudge factor for final descale and range-limit. */
	tmp4 += (RANGE_CENTER << 3) + (1 << 2)

	tmp0 = tmp4 + tmp5
	tmp2 = tmp4 - tmp5

	/* Column 1 */
	tmp4 = DEQUANTIZE(coef_block[DCTSIZE*0+1], qt[quantptr[DCTSIZE*0+1]])
	tmp5 = DEQUANTIZE(coef_block[DCTSIZE*1+1], qt[quantptr[DCTSIZE*1+1]])

	tmp1 = tmp4 + tmp5
	tmp3 = tmp4 - tmp5

	/* Pass 2: process 2 rows, store into output array. */

	/* Row 0 */
	outptr := src[:]

	outptr[0] = range_limit(IRIGHT_SHIFT(tmp0+tmp1, 3))
	outptr[1] = range_limit(IRIGHT_SHIFT(tmp0-tmp1, 3))

	/* Row 1 */
	outptr = src[2:]

	outptr[0] = range_limit(IRIGHT_SHIFT(tmp2+tmp3, 3))
	outptr[1] = range_limit(IRIGHT_SHIFT(tmp2-tmp3, 3))
}

/*
 * Perform dequantization and inverse DCT on one block of coefficients,
 * producing a reduced-size 1x1 output block.
 *
 * We hardly need an inverse DCT routine for this: just take the
 * average pixel value, which is one-eighth of the DC coefficient.
 */

func jpeg_idct_1x1(src, qt *block) {
	var dcval int32
	coef_block := src[:]
	quantptr := jpegZigzagOrder[:]

	/* 1x1 is trivial: just take the DC coefficient divided by 8. */

	dcval = DEQUANTIZE(coef_block[0], qt[quantptr[0]])
	/* Add range center and fudge factor for descale and range-limit. */
	dcval += (RANGE_CENTER << 3) + (1 << 2)

	outptr := src[:]

	outptr[0] = range_limit(IRIGHT_SHIFT(dcval, 3))
}

package jpeg

const CONST_SCALE = ONE << CONST_BITS

/* Convert a positive real constant to an integer scaled by CONST_SCALE.
 * Caution: some C compilers fail to reduce "FIX(constant)" at compile time,
 * thus causing a lot of useless floating-point operations at run time.
 */
func FIX(a float32) int32 {
	return int32(a*CONST_SCALE + 0.5)
}

/*
 * Perform dequantization and inverse DCT on one block of coefficients,
 * producing a reduced-size 7x7 output block.
 *
 * Optimized algorithm with 12 multiplications in the 1-D kernel.
 * cK represents sqrt(2) * cos(K*pi/14).
 */
func idct_7x7(src, qt *block) {
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
	outptr = src[1:]

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

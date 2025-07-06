// Copyright 2011 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package jpegscaled

import (
	"bufio"
	"fmt"
	"image"
	"image/color"
	_ "image/jpeg"
	_ "image/png"
	"os"
	"path/filepath"
	"testing"
)

type imageTest struct {
	goldenFilename string
	filename       string
}

var imageTests = []imageTest{
	// JPEG is a lossy format and hence needs a non-zero tolerance.
	{"testdata/video-001", "testdata/video-001.jpeg"},
	{"testdata/video-001", "testdata/video-001.progressive.jpeg"},
	{"testdata/video-001.221212", "testdata/video-001.221212.jpeg"},
	{"testdata/video-001.cmyk", "testdata/video-001.cmyk.jpeg"},
	{"testdata/video-001.rgb", "testdata/video-001.rgb.jpeg"},
	{"testdata/video-001.progressive.truncated", "testdata/video-001.progressive.truncated.jpeg"},
	{"testdata/video-001.q50.410", "testdata/video-001.q50.410.progressive.jpeg"},
	{"testdata/video-001.q50.411", "testdata/video-001.q50.411.progressive.jpeg"},
	{"testdata/video-001.q50.420", "testdata/video-001.q50.420.progressive.jpeg"},
	{"testdata/video-001.q50.422", "testdata/video-001.q50.422.progressive.jpeg"},
	{"testdata/video-001.q50.440", "testdata/video-001.q50.440.progressive.jpeg"},
	{"testdata/video-001.q50.444", "testdata/video-001.q50.444.progressive.jpeg"},
	{"testdata/video-001.q50.410", "testdata/video-001.q50.410.jpeg"},
	{"testdata/video-001.q50.411", "testdata/video-001.q50.411.jpeg"},
	{"testdata/video-001.q50.420", "testdata/video-001.q50.420.jpeg"},
	{"testdata/video-001.q50.422", "testdata/video-001.q50.422.jpeg"},
	{"testdata/video-001.q50.440", "testdata/video-001.q50.440.jpeg"},
	{"testdata/video-001.q50.444", "testdata/video-001.q50.444.jpeg"},
	{"testdata/video-001.restart2", "testdata/video-001.restart2.jpeg"},
	{"testdata/video-001.separate.dc.progression", "testdata/video-001.separate.dc.progression.progressive.jpeg"},
	{"testdata/video-001.separate.dc.progression", "testdata/video-001.separate.dc.progression.jpeg"},
	// Grayscale images.
	{"testdata/video-005.gray", "testdata/video-005.gray.jpeg"},
	{"testdata/video-005.gray.q50", "testdata/video-005.gray.q50.progressive.jpeg"},
	{"testdata/video-005.gray.q50.2x2", "testdata/video-005.gray.q50.2x2.progressive.jpeg"},
	{"testdata/video-005.gray.q50", "testdata/video-005.gray.q50.jpeg"},
	{"testdata/video-005.gray.q50.2x2", "testdata/video-005.gray.q50.2x2.jpeg"},
}

func TestDecode(t *testing.T) {
	for _, it := range imageTests {
	loop:
		for dctSizeScaled := DCTSIZE; dctSizeScaled > 0; dctSizeScaled-- {
			m, err := decodeJpegScaled(it.filename, dctSizeScaled)
			if err != nil {
				t.Errorf("%s #%d: %v", it.filename, dctSizeScaled, err)
				continue
			}

			goldenFileName := fmt.Sprintf("%s/%s#%d.png", filepath.Dir(it.goldenFilename), filepath.Base(it.goldenFilename),
				dctSizeScaled)

			g, err := decodeStd(goldenFileName)
			if err != nil {
				t.Errorf("decodeStd %s: %v", it.goldenFilename, err)
				continue
			}

			b := g.Bounds()
			if !b.Eq(m.Bounds()) {
				t.Errorf("%s #%d: got bounds %v want %v", it.filename, dctSizeScaled, m.Bounds(), b)
				continue
			}

			for y := b.Min.Y; y < b.Max.Y; y++ {
				for x := b.Min.X; x < b.Max.X; x++ {
					if !withinTolerance(g.At(x, y), m.At(x, y), 1<<8) {
						t.Errorf("%s #%d: at (%d, %d):\ngot  %v\nwant %v",
							it.filename, dctSizeScaled, x, y, rgba(m.At(x, y)), rgba(g.At(x, y)))
						continue loop
					}
				}
			}

			c, _, err := decodeConfig(it.filename)
			if err != nil {
				t.Errorf("%s: %v", it.filename, err)
				continue
			}
			if m.ColorModel() != c.ColorModel {
				t.Errorf("%s #%d: color models differ", it.filename, dctSizeScaled)
				continue
			}
		}
	}
}

func decodeJpegScaled(filename string, dctSizeScaled int) (image.Image, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return Decode(bufio.NewReader(f), dctSizeScaled)
}

func decodeStd(filename string) (image.Image, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	im, _, err := image.Decode(bufio.NewReader(f))
	return im, err
}

func decodeConfig(filename string) (image.Config, string, error) {
	f, err := os.Open(filename)
	if err != nil {
		return image.Config{}, "", err
	}
	defer f.Close()
	return image.DecodeConfig(bufio.NewReader(f))
}

func withinTolerance(c0, c1 color.Color, tolerance int) bool {
	r0, g0, b0, a0 := c0.RGBA()
	r1, g1, b1, a1 := c1.RGBA()
	r := int(delta(r0, r1))
	g := int(delta(g0, g1))
	b := int(delta(b0, b1))
	a := int(delta(a0, a1))
	return r <= tolerance && g <= tolerance && b <= tolerance && a <= tolerance
}

func rgba(c color.Color) string {
	r, g, b, a := c.RGBA()
	return fmt.Sprintf("rgba = 0x%04x, 0x%04x, 0x%04x, 0x%04x for %T%v", r, g, b, a, c, c)
}

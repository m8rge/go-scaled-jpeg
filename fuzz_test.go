// Copyright 2021 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package jpegscaled

import (
	"bytes"
	"fmt"
	"image"
	"image/jpeg"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func FuzzDecode(f *testing.F) {
	if testing.Short() {
		f.Skip("Skipping in short mode")
	}

	testdata, err := os.ReadDir("testdata")
	if err != nil {
		f.Fatalf("failed to read testdata directory: %s", err)
	}
	for _, de := range testdata {
		if de.IsDir() || !strings.HasSuffix(de.Name(), ".jpeg") {
			continue
		}
		b, err := os.ReadFile(filepath.Join("testdata", de.Name()))
		if err != nil {
			f.Fatalf("failed to read testdata: %s", err)
		}
		f.Add(b)
	}

	f.Fuzz(func(t *testing.T, b []byte) {
		cfg, _, err := image.DecodeConfig(bytes.NewReader(b))
		if err != nil {
			return
		}
		if cfg.Width*cfg.Height > 1e6 {
			return
		}
		img, typ, err := image.Decode(bytes.NewReader(b))
		if err != nil || typ != "jpeg" {
			return
		}

		for dctScaledSize := 1; dctScaledSize <= DCTSIZE; dctScaledSize++ {
			t.Run(fmt.Sprintf("dct size %d", dctScaledSize), func(t *testing.T) {
				for q := 1; q <= 100; q++ {
					var w bytes.Buffer
					err := jpeg.Encode(&w, img, &jpeg.Options{Quality: q})
					if err != nil {
						t.Errorf("failed to encode valid image: %s", err)
						continue
					}
					img1, err := Decode(&w, DecodeOptions{DCTSizeScaled: dctScaledSize})
					if err != nil {
						t.Errorf("failed to decode roundtripped image: %s", err)
						continue
					}
					got := img1.Bounds()
					want := img.Bounds()
					want.Max.X *= dctScaledSize / DCTSIZE
					if want.Max.X <= 0 {
						want.Max.X = 1
					}
					want.Max.Y *= dctScaledSize / DCTSIZE
					if want.Max.Y <= 0 {
						want.Max.Y = 1
					}
					if dctScaledSize == DCTSIZE && !got.Eq(want) {
						t.Errorf("image bounds wrong, got: %s, want: %s", got, want)
					}
				}
			})
		}
	})
}

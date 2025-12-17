# go-scaled-jpeg

`go-scaled-jpeg` is a Go library for decoding JPEG images at reduced resolution using scaled IDCT (Inverse Discrete Cosine Transform). It enables memory-efficient image decoding — ideal for thumbnails or previews. Best performance gain achieved with baseline JPEGs. Progressive JPEGs are supported but memory consumption on par with image/jpeg library.

This library combines Go's standard `image/jpeg` infrastructure with a translated and adapted version of the Independent JPEG Group (IJG)'s integer IDCT routines.

## Features

- JPEG decoding at reduced resolutions (1/8, 1/4, 1/2, full)
- Scalable via `DCTSizeScaled` (1–8)
- Tolerant mode to decode incomplete images without failing
- Based on Go standard library and IJG's reference implementation

## Installation

```bash
go get github.com/m8rge/go-scaled-jpeg
```

## Usage

```go
import (
    "os"
    "github.com/m8rge/go-scaled-jpeg"
)

func main() {
    f, err := os.Open("example.jpg")
    if err != nil {
        panic(err)
    }
    defer f.Close()

    img, err := jpegscaled.Decode(f, 4) // Decode at 1/2 resolution
    if err != nil {
        panic(err)
    }

    // Use img as image.Image
}
```

## Benchmarks
```
BenchmarkDecodeBaseline
BenchmarkDecodeBaseline/scaled
BenchmarkDecodeBaseline/scaled/dct_size_1
BenchmarkDecodeBaseline/scaled/dct_size_1-10         	    1905	    544793 ns/op	 907.50 MB/s	   14640 B/op	       5 allocs/op
BenchmarkDecodeBaseline/scaled/dct_size_2
BenchmarkDecodeBaseline/scaled/dct_size_2-10         	    2116	    548515 ns/op	 901.34 MB/s	   16944 B/op	       5 allocs/op
BenchmarkDecodeBaseline/scaled/dct_size_3
BenchmarkDecodeBaseline/scaled/dct_size_3-10         	    2077	    576923 ns/op	 856.96 MB/s	   20656 B/op	       5 allocs/op
BenchmarkDecodeBaseline/scaled/dct_size_4
BenchmarkDecodeBaseline/scaled/dct_size_4-10         	    1988	    595975 ns/op	 829.56 MB/s	   26160 B/op	       5 allocs/op
BenchmarkDecodeBaseline/scaled/dct_size_5
BenchmarkDecodeBaseline/scaled/dct_size_5-10         	    2026	    594736 ns/op	 831.29 MB/s	   32944 B/op	       5 allocs/op
BenchmarkDecodeBaseline/scaled/dct_size_6
BenchmarkDecodeBaseline/scaled/dct_size_6-10         	    1911	    626495 ns/op	 789.15 MB/s	   41136 B/op	       5 allocs/op
BenchmarkDecodeBaseline/scaled/dct_size_7
BenchmarkDecodeBaseline/scaled/dct_size_7-10         	    1825	    660812 ns/op	 748.17 MB/s	   54832 B/op	       5 allocs/op
BenchmarkDecodeBaseline/scaled/dct_size_8
BenchmarkDecodeBaseline/scaled/dct_size_8-10         	    1664	    705304 ns/op	 700.97 MB/s	   63024 B/op	       5 allocs/op
BenchmarkDecodeBaseline/optimized
BenchmarkDecodeBaseline/optimized-10                 	    1684	    721496 ns/op	 685.24 MB/s	   63024 B/op	       5 allocs/op
BenchmarkDecodeProgressive
BenchmarkDecodeProgressive/scaled
BenchmarkDecodeProgressive/scaled/dct_size_1
BenchmarkDecodeProgressive/scaled/dct_size_1-10         	    1320	    875977 ns/op	 564.40 MB/s	  211337 B/op	      13 allocs/op
BenchmarkDecodeProgressive/scaled/dct_size_2
BenchmarkDecodeProgressive/scaled/dct_size_2-10         	    1300	    873480 ns/op	 566.01 MB/s	  213632 B/op	      13 allocs/op
BenchmarkDecodeProgressive/scaled/dct_size_3
BenchmarkDecodeProgressive/scaled/dct_size_3-10         	    1383	    879395 ns/op	 562.20 MB/s	  217344 B/op	      13 allocs/op
BenchmarkDecodeProgressive/scaled/dct_size_4
BenchmarkDecodeProgressive/scaled/dct_size_4-10         	    1347	    891384 ns/op	 554.64 MB/s	  222848 B/op	      13 allocs/op
BenchmarkDecodeProgressive/scaled/dct_size_5
BenchmarkDecodeProgressive/scaled/dct_size_5-10         	    1326	    913840 ns/op	 541.01 MB/s	  229632 B/op	      13 allocs/op
BenchmarkDecodeProgressive/scaled/dct_size_6
BenchmarkDecodeProgressive/scaled/dct_size_6-10         	    1276	    922707 ns/op	 535.81 MB/s	  237824 B/op	      13 allocs/op
BenchmarkDecodeProgressive/scaled/dct_size_7
BenchmarkDecodeProgressive/scaled/dct_size_7-10         	    1239	    961199 ns/op	 514.36 MB/s	  251520 B/op	      13 allocs/op
BenchmarkDecodeProgressive/scaled/dct_size_8
BenchmarkDecodeProgressive/scaled/dct_size_8-10         	    1156	   1028548 ns/op	 480.68 MB/s	  259712 B/op	      13 allocs/op
BenchmarkDecodeProgressive/optimized
BenchmarkDecodeProgressive/optimized-10                 	    1147	   1043039 ns/op	 474.00 MB/s	  259712 B/op	      13 allocs/op
```

## License

This project contains code under multiple licenses:

- MIT License: For original code written for this project (see LICENSE).
- IJG License: For translated parts from the Independent JPEG Group (e.g., IDCT routines) (see COPYING.IJG).
- BSD 3-Clause License: For parts derived from the Go standard library (see LICENSE.GOLANG).

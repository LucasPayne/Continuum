#!/bin/bash

convert +append quadratic* tmp_1.png
convert +append linear* tmp_2.png
convert -append tmp_1.png tmp_2.png collage_1.png


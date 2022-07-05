#!/bin/bash
# verschiebt die von Latex und biber erzeugten Zwischendateien nach /tmp
mv *.aux *.bbl *.bcf *.blg *.log *.out *.run.xml *.toc *.synctex.gz /tmp

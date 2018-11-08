#from chen
nucmer -c 200 -g 200 -p out allpathslg.fasta soapdenovo.fasta
show-coords -c -d -l -I 95 -L 10000 -r out.delta > out.show
mummerplot -f -l -p out -s large -t png out.delta
mummerplot -f -l -p out1 -s large -t png -r scaffold_0 out.delta
mummerplot -f -l -p out2 -s large -t png -S out.delta
/opt/pacbio_soft/smrtanalysis/userdata/jobs/016/016938

/opt/pacbio_soft/smrtanalysis/current/smrtcmds/bin

smrtpipe.py -D NPROC=3 -D CLUSTER=BASH -D MAX_THREADS=4 --params=$TESTDIR/params.xml xml:$TESTDIR/input.xml >smrtpipe.log

TESTDIR=/opt/pacbio_soft/smrtanalysis/current/doc/examples/smrtpipe_assembly_hgap3

/opt/pacbio_soft/smrtanalysis/current/analysis/bin/fofnToSmrtpipeInput.py input.fofn > input.xml

#用smrtlink的python 2.7.9
nohup /opt/pacbio_soft/smrtanalysis/current/smrtcmds/bin/smrtpipe -D NPROC=3 -D CLUSTER=BASH -D MAX_THREADS=20 --params=/opt/pacbio_soft/smrtanalysis/current/doc/examples/smrtpipe_assembly_hgap3/params.xml xml:input.xml --output=./ >smrtpipe.log &


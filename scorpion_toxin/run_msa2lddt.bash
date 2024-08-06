#!/bin/bash -e

msa=aln.fa
files=pdbfiles.txt

if [ x$files == x ] ; then
	echo Missing arg >> /dev/stderr
	exit 1
fi

if [ ! -s $msa ] ; then
	echo Not found msa=$msa >> /dev/stderr
fi

if [ ! -s $files ] ; then
	echo Not found files=$files >> /dev/stderr
fi

foldmason=/z/sw/foldmason/foldmason/bin/foldmason

tmpdir=foldmason_tmp
pdbdir=foldmason_pdbdir

rm -rf $tmpdir $pdbdir
mkdir -p $tmpdir $pdbdir

cp -v `cat $files` $pdbdir
cp -v $msa $tmpdir/msa.fa

$foldmason createdb $pdbdir tmpdb
$foldmason msa2lddt tmpdb $tmpdir/msa.fa

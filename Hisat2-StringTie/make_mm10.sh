#!/bin/sh

#
# Downloads sequence for the MM10 version of M. musculus (mouse) from
# UCSC.
#
# The base files, named ??.fa.gz
#
# By default, this script builds and index for just the base files,
# since alignments to those sequences are the most useful.  To change
# which categories are built by this script, edit the CHRS_TO_INDEX
# variable below.
#

UCSC_MM10_BASE=http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips
F=chromFa.tar.gz

get() {
	file=$1
	if ! wget --version >/dev/null 2>/dev/null ; then
		if ! curl --version >/dev/null 2>/dev/null ; then
			echo "Please install wget or curl somewhere in your PATH"
			exit 1
		fi
		curl -o `basename $1` $1
		return $?
	else
		wget $1
		return $?
	fi
}

HISAT2_BUILD_EXE=./hisat2-build
if [ ! -x "$HISAT2_BUILD_EXE" ] ; then
	if ! which hisat2-build ; then
		echo "Could not find hisat2-build in current directory or in PATH"
		exit 1
	else
		HISAT2_BUILD_EXE=`which hisat2-build`
	fi
fi

rm -f genome.fa
get ${UCSC_MM10_BASE}/$F || (echo "Error getting $F" && exit 1)
tar xvzfO $F > genome.fa || (echo "Error unzipping $F" && exit 1)
rm $F

CMD="${HISAT2_BUILD_EXE} genome.fa genome"
echo Running $CMD
if $CMD ; then
	echo "genome index built; you may remove fasta files"
else
	echo "Index building failed; see error message"
fi

/bin/bash /app/readme.hisat2.stringtie.DESEQ2.txt
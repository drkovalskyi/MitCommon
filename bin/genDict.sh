#!/bin/bash

if ! [ $CMSSW_BASE ] || ! [ $SCRAM_ARCH ]
then
  echo "CMSSW_BASE and SCRAM_ARCH must be set."
  exit 1
fi

PACKAGES="$@"

echo "Generating ROOT dictionaries for:"
echo " $PACKAGES"

for PACKAGE in $PACKAGES
do
  DICTDIR=$CMSSW_BASE/src/$PACKAGE/dict
  SRCDIR=$CMSSW_BASE/src/$PACKAGE/src
  INCDIR=$CMSSW_BASE/src/$PACKAGE/interface

  LIBDIR=$CMSSW_BASE/lib/$SCRAM_ARCH
  [ -d $LIBDIR ] || mkdir $LIBDIR

  cd $SRCDIR

  for LINKDEF in $(ls $DICTDIR/*LinkDef.h); do
    # collect include headers from LinkDef file
    HEADERS=$(sed -n 's|#include *"\([^"]*\)"|'$CMSSW_BASE'/src/\1|p' $LINKDEF | tr '\n' ' ')
    # the name of the dictionary code file generated from ABCLinkDef.h will be ABC_LinkDefDict.cc
    OUTPUT=$(sed 's/LinkDef.h/_LinkDefDict/' <<< $(basename $LINKDEF))

    echo rootcling -f $OUTPUT.cc -I$CMSSW_BASE/src $HEADERS $LINKDEF
    rootcling -f $OUTPUT.cc -I$CMSSW_BASE/src $HEADERS $LINKDEF

    echo mv ${OUTPUT}_rdict.pcm $LIBDIR/
    mv ${OUTPUT}_rdict.pcm $LIBDIR/
  done
  
  if [ -e $DICTDIR/classes_def.xml ]
  then
    # the name of the dictionary code file generated from classes_xml will be PACKAGE_ReflexDict.cc
    OUTPUT=$(sed 's|\([^/]*\)/\(.*\)|\1\2|' <<< $PACKAGE)_ReflexDict

    echo genreflex $DICTDIR/classes.h -s $DICTDIR/classes_def.xml -o $OUTPUT.cc --rootmap=$OUTPUT.rootmap -I$CMSSW_BASE/src
    genreflex $DICTDIR/classes.h -s $DICTDIR/classes_def.xml -o $OUTPUT.cc --rootmap=$OUTPUT.rootmap -I$CMSSW_BASE/src

    echo mv $OUTPUT.rootmap $LIBDIR/
    mv $OUTPUT.rootmap $LIBDIR/
    echo mv ${OUTPUT}_rdict.pcm $LIBDIR/
    mv ${OUTPUT}_rdict.pcm $LIBDIR/
  fi
done

echo "=== genDict.sh ==="
echo "ROOT dictionary binaries were copied to $CMSSW_BASE/lib/$SCRAM_ARCH."
echo "Do not forget to run"
echo '$CMSSW_BASE/src/MitCommon/bin/genDict.sh [packages]'
echo "after next scram build clean."
echo "=================="
echo ""

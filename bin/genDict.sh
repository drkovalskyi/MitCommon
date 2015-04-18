#!/bin/bash

check-update () {
  local DEPENDANT=$1
  local SOURCES="$@"

  local UPDATED=false
    
  if ! [ -e $DEPENDANT ]
  then
    UPDATED=true
  else
    LASTUPDATE=$(stat -c %Z $DEPENDANT)
    local SOURCE
      
    for SOURCE in $SOURCES
    do
      if [ $LASTUPDATE -lt $(stat -c %Z $SOURCE) ]
      then
        UPDATED=true
        break
      fi
    done
  fi

  echo $UPDATED
}


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
    for HEADER in $HEADERS
    do
      if ! [ -e $HEADER ]
      then
        echo "$HEADER does not exist"
        exit 1
      fi
    done
      
    # the name of the dictionary code file generated from ABCLinkDef.h will be ABC_LinkDefDict.cc
    OUTPUT=$(sed 's/LinkDef.h/_LinkDefDict/' <<< $(basename $LINKDEF))

    # generate dictionary if $OUTPUT.cc does not exist or is older than one of the source files

    if $(check-update $OUTPUT.cc $HEADERS $LINKDEF)
    then
      echo rootcling -f $OUTPUT.cc -I$CMSSW_BASE/src $HEADERS $LINKDEF
      rootcling -f $OUTPUT.cc -I$CMSSW_BASE/src $HEADERS $LINKDEF

      echo mv ${OUTPUT}_rdict.pcm $LIBDIR/
      mv ${OUTPUT}_rdict.pcm $LIBDIR/
    fi
  done
  
  if [ -e $DICTDIR/classes_def.xml ]
  then
    # the name of the dictionary code file generated from classes_xml will be PACKAGE_ReflexDict.cc
    OUTPUT=$(sed 's|\([^/]*\)/\(.*\)|\1\2|' <<< $PACKAGE)_ReflexDict

    if $(check-update $OUTPUT.cc $DICTDIR/classes_def.xml $DICTDIR/classes.h)
    then
      echo genreflex $DICTDIR/classes.h -s $DICTDIR/classes_def.xml -o $OUTPUT.cc --rootmap=$OUTPUT.rootmap -I$CMSSW_BASE/src
      genreflex $DICTDIR/classes.h -s $DICTDIR/classes_def.xml -o $OUTPUT.cc --rootmap=$OUTPUT.rootmap -I$CMSSW_BASE/src

      echo mv $OUTPUT.rootmap $LIBDIR/
      mv $OUTPUT.rootmap $LIBDIR/
      echo mv ${OUTPUT}_rdict.pcm $LIBDIR/
      mv ${OUTPUT}_rdict.pcm $LIBDIR/
    fi
  fi
done

echo "=== genDict.sh ==="
echo "ROOT dictionary binaries were copied to $CMSSW_BASE/lib/$SCRAM_ARCH."
echo "Do not forget to run"
echo '$CMSSW_BASE/src/MitCommon/bin/genDict.sh [packages]'
echo "after next scram build clean."
echo "=================="
echo ""

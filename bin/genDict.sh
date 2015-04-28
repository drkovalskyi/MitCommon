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
      if ! [ -e $SOURCE ]
      then
        echo "$SOURCE does not exist"
        exit 1
      fi
      
      if [ $LASTUPDATE -lt $(stat -c %Z $SOURCE) ]
      then
        UPDATED=true
        break
      fi
    done
  fi

  echo $UPDATED
}

usage () {
  echo "Usage: genDict.sh [-c] [-f] PACKAGES"
  echo "Generate ROOT dictionary source code. PACKAGE is e.g. MitAna/DataTree."
  echo ""
  echo "  -c  Clear generated source code."
  echo "  -f  Force generation. Otherwise dictionaries are written only if"
  echo "      header files are updated."

  exit $1
}


FORCE=false
CLEAR=false
while getopts fch OPT; do
  case $OPT in
    c)
      CLEAR=false
      ;;
        
    f)
      FORCE=true
      ;;
    h)
      usage 0
      ;;
    \?)
      echo " Invalid option: -$OPTARG" >& 2
      usage 1
      ;;
  esac
done

shift $((OPTIND-1))

if ! [ $CMSSW_BASE ] || ! [ $SCRAM_ARCH ]
then
  echo "CMSSW_BASE and SCRAM_ARCH must be set."
  exit 1
fi

PACKAGES="$@"

if ! [ "$PACKAGES" ]
then
  for DIR in $(ls $CMSSW_BASE/src)
  do
    for SUBDIR in $(ls $CMSSW_BASE/src/$DIR)
    do
      [ -d $CMSSW_BASE/src/$DIR/$SUBDIR/dict ] && PACKAGES="$PACKAGES $DIR/$SUBDIR"
    done
  done
fi

echo "Generating ROOT dictionaries for:"
echo " $PACKAGES"

for PACKAGE in $PACKAGES
do
  DICTDIR=$CMSSW_BASE/src/$PACKAGE/dict
  SRCDIR=$CMSSW_BASE/src/$PACKAGE/src
  INCDIR=$CMSSW_BASE/src/$PACKAGE/interface

  if ! [ -d $DICTDIR ] || ! [ -d $SRCDIR ] || ! [ -d $INCDIR ]
  then
    echo "$PACKAGE does not appear to be a valid package."
    continue
  fi

  LIBDIR=$CMSSW_BASE/lib/$SCRAM_ARCH
  [ -d $LIBDIR ] || mkdir $LIBDIR

  MAKEFILE=$DICTDIR/Makefile

  cd $SRCDIR

  for LINKDEF in $(ls $DICTDIR/*LinkDef.h); do
    # the name of the dictionary code file generated from ABCLinkDef.h will be ABC_LinkDefDict.cc
    OUTPUT=$(sed 's/LinkDef.h/_LinkDefDict/' <<< $(basename $LINKDEF))

    if $CLEAR
    then
      rm -f $OUTPUT 2>/dev/null
      rm -f $LIBDIR/${OUTPUT}_rdict.pcm 2>/dev/null
    else
      INCLUDES=$(makedepend -I$CMSSW_BASE/src -f- $LINKDEF 2>/dev/null | awk '/Mit/ {print $2}' | tr '\n' ' ')
        
      # generate dictionary if $OUTPUT.cc does not exist or is older than one of the source files
      if $FORCE || $(check-update $OUTPUT.cc $INCLUDES $LINKDEF) || ! [ -e $LIBDIR/${OUTPUT}_rdict.pcm ]
      then
        HEADERS=$(sed -n 's|#include *"\([^"]*\)"|'$CMSSW_BASE'/src/\1|p' $LINKDEF | tr '\n' ' ')
          
        echo rootcling -f $OUTPUT.cc -I$CMSSW_BASE/src $HEADERS $LINKDEF
        rootcling -f $OUTPUT.cc -I$CMSSW_BASE/src $HEADERS $LINKDEF
  
        echo mv ${OUTPUT}_rdict.pcm $LIBDIR/
        mv ${OUTPUT}_rdict.pcm $LIBDIR/
      fi
    fi
  done
  
  if [ -e $DICTDIR/classes_def.xml ]
  then
    # the name of the dictionary code file generated from classes_xml will be PACKAGE_ReflexDict.cc
    OUTPUT=$(sed 's|\([^/]*\)/\(.*\)|\1\2|' <<< $PACKAGE)_ReflexDict

    if $CLEAR
    then
      rm -f $OUTPUT 2>/dev/null
      rm -f $LIBDIR/$OUTPUT.rootmap 2>/dev/null
      rm -f $LIBDIR/${OUTPUT}_rdict.pcm 2>/dev/null
    else
      INCLUDES=$(makedepend -I$CMSSW_BASE/src -f- $DICTDIR/classes.h 2>/dev/null | awk '/Mit/ {print $2}' | tr '\n' ' ')      

      if $FORCE || $(check-update $OUTPUT.cc $DICTDIR/classes_def.xml $DICTDIR/classes.h $INCLUDES) ||
              ! [ -e $LIBDIR/$OUTPUT.rootmap ] || ! [ -e $LIBDIR/${OUTPUT}_rdict.pcm ]
      then
        echo genreflex $DICTDIR/classes.h -s $DICTDIR/classes_def.xml -o $OUTPUT.cc --rootmap=$OUTPUT.rootmap -I$CMSSW_BASE/src
        genreflex $DICTDIR/classes.h -s $DICTDIR/classes_def.xml -o $OUTPUT.cc --rootmap=$OUTPUT.rootmap -I$CMSSW_BASE/src
  
        echo mv $OUTPUT.rootmap $LIBDIR/
        mv $OUTPUT.rootmap $LIBDIR/
        echo mv ${OUTPUT}_rdict.pcm $LIBDIR/
        mv ${OUTPUT}_rdict.pcm $LIBDIR/
      fi
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

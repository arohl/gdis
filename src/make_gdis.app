#!/bin/bash

################################################################################
# Environment variable check.  These should be set in the makefile.

if [[ -z $PLATYPUS ]]; then
  echo "Error: environment variable PLATYPUS does not exist"
  echo "I hope you are running this script from the makefile."
  exit 42
fi
if [[ -z $FINK ]]; then
  echo "Error: environment variable FINK does not exist"
  echo "I hope you are running this script from the makefile."
  exit 42
fi
if [[ ! -d $FINK ]]; then
  echo "Error: environment variable FINK does not point to a real directory."
  echo "$FINK"
  exit 42
fi

################################################################################
# A convenience function.

function checkfile() {
  if [[ ! -f $1 ]]; then
    echo "Error: file $1 does not exist"
    exit 69
  fi
}

################################################################################

# $FINK with slashes escaped.
FINK1=$( echo $FINK | sed 's/\//\\\//g' )

rm -rf gdis.app

checkfile gdis
checkfile gdis.elements
checkfile gdis.manual
checkfile gdis.library
checkfile GDIS.icns
checkfile gdis.app.template/Contents/Resources/script

echo "Running Platypus..."
$PLATYPUS -a GDIS -f gdis -f gdis.elements -f gdis.manual -f gdis.library -t python -o None -i GDIS.icns -V 0.89 -u "Sean Fleming" -I "org.sean.gdis" -D gdis.app.template/Contents/Resources/script gdis.app
if [[ ! -f gdis.app/Contents/Resources/gdis ]]; then
  echo "Error: Platypus failed."
  exit 1
else 
  echo "done."
fi

# Fix in case script is not executable after above step.
if [[ ! -x gdis.app/Contents/Resources/script ]]; then
  chmod a+x gdis.app/Contents/Resources/script
fi

# Workaround for platypus icon copy bug
if [[ ! -f gdis.app/Contents/Resources/GDIS.icns ]]; then
  cp GDIS.icns gdis.app/Contents/Resources
  sed -i "" 's/appIcon.icns/GDIS.icns/' gdis.app/Contents/Info.plist
fi

mkdir -p gdis.app/Contents/Frameworks
#Make fink library dependencies to be within gdis.app.
echo "Fixing library dependencies in gdis..."
for lib in `otool -L gdis.app/Contents/Resources/gdis | grep -i "$FINK"` ; do
  if [[ -n $( echo $lib | grep -i $FINK ) ]]; then
    checkfile $lib
    cp -f $lib gdis.app/Contents/Frameworks
    newlib=$( echo "$lib" | sed 's/^.*\//\@executable_path\/..\/Frameworks\//' )
    echo "   $newlib"
    install_name_tool -change $lib $newlib gdis.app/Contents/Resources/gdis
  fi
done

echo "Copying pango..."
for file in `find $FINK/lib/pango` ; do
#  newfile=$( echo "$file" | sed "s/$FINK1\/lib/gdis.app\/Contents\/Frameworks/" )
  newfile=$( echo "$file" | perl -p -e "s/$FINK1\/lib/gdis.app\/Contents\/Frameworks/gi" )
  if [[ -d $file ]]; then
    mkdir -p $newfile
  else
    newfile1=$( echo "$newfile" | grep -E -v -e "\.(la|a)$" )
    if [[ -n "$newfile1" ]] ; then
      echo "  $newfile1"
      checkfile $file
      cp $file $newfile
    fi
  fi
done

#Some fink libraries are dependent on other fink libraries.  For now copy them
#manually.
echo "Copying extra libs to gdis.app/Contents/Frameworks"
checkfile $FINK/lib/libintl.3.dylib
cp $FINK/lib/libintl.3.dylib gdis.app/Contents/Frameworks

#Got to copy pango.modules to Resources/etc/pango and put in $CWD/lib/pango/...
#Also copy pango.aliases there too.
mkdir -p gdis.app/Contents/Resources/etc/pango
checkfile $FINK/etc/pango/pangox.aliases
cp $FINK/etc/pango/pangox.aliases gdis.app/Contents/Resources/etc/pango/pangox.aliases
#pango-querymodules | grep -v "^#" | sed "s/$FINK1\/lib/\$\{CWD\}\/..\/Frameworks/" | sed -E "s/^([^ ]+)/\"\1\"/" > gdis.app/Contents/Resources/etc/pango/pango.modules
pango-querymodules | grep -v "^#" | perl -p -e "s/$FINK1\/lib/\\$\{\CWD}\/..\/Frameworks/gi" | sed -E "s/^([^ ]+)/\"\1\"/" > gdis.app/Contents/Resources/etc/pango/pango.modules

echo "Copying GTK2..."
mkdir -p gdis.app/Contents/Frameworks/gtk-2.0/$GTK2VERS
for file in `find $FINK/lib/gtk-2.0/$GTK2VERS` ; do
#  newfile=$( echo "$file" | sed "s/$FINK1\/lib/gdis.app\/Contents\/Frameworks/" )
  newfile=$( echo "$file" | perl -p -e "s/$FINK1\/lib/gdis.app\/Contents\/Frameworks/gi" )
  if [[ -d $file ]]; then
    mkdir -p $newfile
  else
    newfile1=$( echo "$newfile" | grep -E -v -e "\.(la|a)$" )
    if [[ -n "$newfile1" ]] ; then
      echo "   $newfile1"
      checkfile $file
      cp $file $newfile
    fi
  fi
done

mkdir -p gdis.app/Contents/Resources/etc/gtk-2.0
#gtk-query-immodules-2.0 | sed "s/$FINK1\/lib/\$\{CWD\}\/..\/Frameworks/" | grep -v "^#" > gdis.app/Contents/Resources/etc/gtk-2.0/gtk.immodules
gtk-query-immodules-2.0 | perl -p -e "s/$FINK1\/lib/\\$\{\CWD}\/..\/Frameworks/gi" | grep -v "^#" > gdis.app/Contents/Resources/etc/gtk-2.0/gtk.immodules

# Get rid of fink dependencies out of libraries too.
echo "Fixing inter-library dependencies..."
for file in `find gdis.app/Contents/Frameworks -name "*.dylib" -or -name "*.so"`; do
  echo "  processing $file"
  newlib=$( echo "$file" | sed 's/^.*\//\@executable_path\/..\/Frameworks\//' )
  install_name_tool -id $newlib $file
  for lib in `otool -L $file | grep -i "$FINK"` ; do
    if [[ -n $( echo $lib | grep -i $FINK ) ]]; then
      newlib=$( echo "$lib" | sed 's/^.*\//\@executable_path\/..\/Frameworks\//' )
      install_name_tool -change $lib $newlib $file
    fi
  done
done

gdk-pixbuf-query-loaders | perl -p -e "s/$FINK1\/lib/\\$\{\CWD}\/..\/Frameworks/gi"> gdis.app/Contents/Resources/gdk-pixbuf.loaders
#gdk-pixbuf-query-loaders | sed "s/$FINK1\/lib/\$\{CWD\}\/..\/Frameworks/"> gdis.app/Contents/Resources/gdk-pixbuf.loaders

# Final check that we've sorted out all the Fink dependencies.
echo "Just checking sure we got all the Fink dependencies sorted..."
for file in `find gdis.app/Contents/Frameworks -name "*.dylib" -or -name "*.so"`; do
  if [[ -n `otool -L $file | grep -i "$FINK"` ]]; then
    echo "$file still has a FINK depenency"
  fi
done
if [[ -n `otool -L gdis.app/Contents/Resources/gdis | grep -i "$FINK"` ]]; then
  echo "gdis.app/Contents/Resources/gdis still has a FINK depenency"
fi


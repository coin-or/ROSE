#!/bin/sh

set -e

wgetcmd=wget
wgetcount=`which wget 2>/dev/null | wc -w`
if test ! $wgetcount = 1; then
  echo "Utility wget not found in your PATH."
  if test `uname` = Darwin; then
    wgetcmd=ftp
    echo "Using ftp command instead."
  else
    exit -1
  fi
fi

echo " "
echo "Running script for downloading the source code for the ASL"
echo " "

if [ ! -d amplsolver ] ; then

  if [ -f solvers.tar.gz ] ; then
    gunzip solvers.tar.gz 
  else 
    echo "Downloading the source code from www.netlib.org..."
    $wgetcmd ftp://www.netlib.org/ampl/solvers.tar
  fi
  rm -rf solvers

  echo "Unpacking the source code..."
  tar xf solvers.tar

  echo "Uncompressing the source files..."
  gunzip -fr solvers

  echo "Adding No_dtoa to CFLAGS..."
  cd solvers
  mv makefile.u makefile.u.orig
  sed -e 's/CFLAGS = /CFLAGS = -DNo_dtoa /g' makefile.u.orig > makefile.u 
  cd ..

  echo "Compressing the tar file..."
  gzip solvers.tar
  mv solvers amplsolver
fi

cd amplsolver/
sed -e 's/`uname -i` 2>\/dev\/null/`uname -i 2>\/dev\/null`/g' configurehere > ttt
mv ttt configurehere
chmod +x configure
chmod +x configurehere
./configurehere
sed -i "s/^CFLAGS\ =\ -DNo_dtoa\ -O$/CFLAGS\ =\ -O/" makefile

echo " "
echo "Done downloading the source code for ASL."
echo " "


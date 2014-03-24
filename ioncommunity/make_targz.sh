#!/bin/bash

PLUGIN_DIR=AssemblerPlus
rm -rf $PLUGIN_DIR/bin/mira* $PLUGIN_DIR/bin/SPAdes* $PLUGIN_DIR/bin/quast*
mkdir -p $PLUGIN_DIR/bin

#SPADES_REPO_DIR=$HOME/algorithmic-biology
QUAST_REPO_DIR=$HOME/quast

# For building a plugin with different versions of software:
# 1) edit select options in instance.html
# 2) bump version in launch.sh
# 3) when QUAST is updated, edit $quastPath in bin/assembler.pl

echo "Downloading and extracting MIRA 4.0..."
wget -q http://sourceforge.net/projects/mira-assembler/files/MIRA/stable/mira_4.0_linux-gnu_x86_64_static.tar.bz2/download && tar xjf download && rm download
mv mira_4.0_linux-gnu_x86_64_static/ $PLUGIN_DIR/bin/mira-4.0

for ver in "2.5.1" "3.0.0"; do
echo "Downloading and extracting SPAdes $ver..."
wget -q http://spades.bioinf.spbau.ru/release${ver}/SPAdes-${ver}-Linux.tar.gz
tar xzf SPAdes-${ver}-Linux.tar.gz
mv SPAdes-${ver}-Linux $PLUGIN_DIR/bin
rm SPAdes-${ver}-Linux.tar.gz
done

echo "Creating Quast archive using make_tar.sh..."
cd $QUAST_REPO_DIR && bash make_tar.sh && cd -
QUAST_VERSION=`cat $QUAST_REPO_DIR/VERSION`
tar xzf $QUAST_REPO_DIR/quast-$QUAST_VERSION.tar.gz 
mv quast-$QUAST_VERSION/ $PLUGIN_DIR/bin
echo "Copied Quast files to $PLUGIN_DIR/bin/ directory"

# download DownsampleSam.jar from an arbitrary place (to avoid downloading extra ~50MB of other tools from SourceForge)
wget -q --no-check-certificate http://github.com/fgvieira/ngsClean/raw/ead0d60af81edcbbcef4e655d6a44aa360a1a6cc/external_progs/picard-tools-1.90/DownsampleSam.jar
mv DownsampleSam.jar $PLUGIN_DIR/bin
echo "Downloaded bin/DownsampleSam.jar to $PLUGIN_DIR/bin/"

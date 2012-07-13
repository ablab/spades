sudo apt-get install build-essential bzip2 libbz2-dev libz-dev

ftp_proxy="http://192.168.0.2:3128" wget ftp://ftp.gnu.org/gnu/gmp/gmp-5.0.5.tar.bz2
wget http://mirrors-us.seosue.com/gcc/releases/gcc-4.4.7/gcc-4.4.7.tar.gz
wget http://www.mpfr.org/mpfr-current/mpfr-3.1.1.tar.gz
wget http://www.multiprecision.org/mpc/download/mpc-0.9.tar.gz
wget http://www.bastoul.net/cloog/pages/download/count.php3?url=./cloog-0.17.0.tar.gz
wget http://bugseng.com/products/ppl/download/ftp/releases/0.12.1/ppl-0.12.1.tar.gz
wget http://downloads.sourceforge.net/project/boost/boost/1.49.0/boost_1_49_0.tar.gz
wget http://www.cmake.org/files/v2.6/cmake-2.6.4.tar.gz
bzip2 -d gmp-5.0.5.tar.bz2
tar -xf gmp-5.0.5.tar
tar -xzf mpfr-3.1.1.tar.gz
tar -xzf gcc-4.4.7.tar.gz
tar -xzf mpc-0.9.tar.gz
tar -xzf cloog-0.17.0.tar.gz
tar -xzf ppl-0.12.1.tar.gz
tar -xzf boost_1_49_0.tar.gz
tar -xzf cmake-2.6.4.tar.gz


cd gmp-5.0.5/
./configure --enable-cxx
make
sudo make install

cd ../mpfr-3.1.0
./configure
make
sudo make install

cd ../mpc-0.9
./configure
make
sudo make install

cd ../cloog-0.17.0
./configure
make
sudo make install

cd ../ppl-0.12.1/
./configure
make
sudo make install

cd ../gcc-4.4.7/
export LD_LIBRARY_PATH=/usr/local/lib/
./configure --enable-languages=c,c++ --disable-multilib
./configure --enable-languages=c,c++ --disable-multilib --disable-ppl-version-check  --disable-cloog-version-check
make
sudo make install

cd ../boost_1_49_0
./bootstrap.sh --with-libraries=system,filesystem,iostreams,serialization,thread
./b2
sudo ./b2 install

cd ../cmake-2.6.4/
CC=gcc CXX=g++ ./bootstrap
make
sudo make install

cd ..
sudo apt-get install tcl8.4 gettext libcurl3-openssl-dev


wget http://git-core.googlecode.com/files/git-1.7.11.1.tar.gz
tar -xzf git-1.7.11.1.tar.gz
cd git-1.7.11.1/
./configure --prefix=/usr
make
sudo make install

env GIT_SSL_NO_VERIFY=true git clone https://github.com/ablab/algorithmic-biology.git
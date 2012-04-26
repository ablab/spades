rm -r rpmbuild
rpmdev-setuptree
scp builder@192.168.222.254:release2.0.1/spades.spec .
cp spades.spec rpmbuild/SPECS
cd rpmbuild/SOURCES   
scp builder@192.168.222.254:release2.0.1/spades_2.0.1.tar.gz .
tar -xvf spades_2.0.1.tar.gz
cd spades-2.0.1
touch configure
chmod a+x configure 
cd ..
rm spades_2.0.1.tar.gz
tar -czvf spades_2.0.1.tar.gz spades-2.0.1
cd ../SPECS
rpmbuild -ba spades.spec
scp ../RPMS/x86_64/spades-2.0.1-1.fc16.x86_64.rpm builder@192.168.222.254:rpm
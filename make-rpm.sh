 rpmdev-setuptree
 cp spades.spec rpmbuild/SPECS
 cd rpmbuild/SOURCES   
 scp yasha@192.168.222.223:spades2/spades_2.0.0.tar.gz .
 tar -xvf spades_2.0.0.tar.gz
 cd spades-2.0.0
 touch configure
 chmod a+x configure 
 cd ..
 rm spades_2.0.0.tar.gz
 tar -czvf spades_2.0.0.tar.gz spades-2.0.0
 cd ../SPECS
 rpmbuild -ba spades.spec
 scp ../RPMS/x86_64/spades-2.0.0-1.fc16.x86_64.rpm yasha@192.168.222.223:rpm
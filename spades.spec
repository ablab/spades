%define debug_package %{nil}
Name:           spades 
Version:        %version
Release:        1%{?dist}
License:        GPLv2
Summary:        SPAdes genome assembler
URL:            http://bioinf.spbau.ru/spades
BuildRoot:      %{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)
Source: 	spades_2.0.1.tar.gz
Requires:       cmake boost-devel python-devel zlib-devel
Autoreq:        0
%description  

%prep
%setup -q

%build

%configure

%install
rm -rf %{buildroot}
make install DESTDIR=%{buildroot}

%clean
rm -rf %{buildroot}

%files
%defattr(-,root,root,-)
/usr/share/spades
/usr/bin/spades*.py
%doc

%changelog
FROM ubuntu:20.04

MAINTAINER GlaDIAtorAdmin

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y apt-utils

RUN apt-get update \
    && apt-get -y install wget\
    && apt-get -y install git\
    && apt-get -y install emacs-nox\
    && apt-get -y install build-essential\
    && apt-get -y install openssh-server\
    && apt-get -y install sshfs\
    && apt-get -y install tandem-mass\
    && apt-get -y install openjdk-8-jdk\
    && apt-get -y install screen\
    && apt-get -y install xpra\
    && apt-get -y install g++\
    && apt-get -y install zlib1g-dev\
    && apt-get -y install libghc-bzlib-dev\
    && apt-get -y install gnuplot\
    && apt-get -y install unzip\
    && apt-get -y install locales\
    && apt-get -y install expat\
    && apt-get -y install libexpat1-dev\
    && apt-get -y install subversion\
    && apt-get -y install comet-ms\
    && apt-get -y install libfindbin-libs-perl\
    && apt-get -y install libxml-parser-perl\
    && apt-get -y install libtool-bin\
    && apt-get -y install curl\
    && apt-get -y install sudo \
    && apt-get -y install cmake\
    && apt-get -y install gfortran-multilib\
    && apt-get -y install python3-plotly\
    && apt-get -y install python3-pandas \
    && apt-get -y install python3-pandas-lib \
    && apt-get -y install python3-pip\
    && apt-get -y install python3-pymzml\
    && apt-get -y install python3-psutil\
    && apt-get -y install python3-virtualenv\
    && apt-get -y install python3-pyramid\
    && apt-get -y install python-numpy\
    && apt-get -y install npm \
    && apt-get -y install mono-complete \
    && apt-get -y install python3-biopython \
    && apt-get -y install libpwiz3 \
    && apt-get -y install libpwiz-dev \
    && apt-get -y install libpwiz-tools \
    && apt-get -y install openms \
    && apt-get -y install libgd-dev \
    && apt-get -y install cython3 \
    && apt-get -y install python-dev \
    && apt-get -y install libxml2-dev \
    && apt-get -y install libcurl4-openssl-dev


RUN apt-get -y purge r-base* r-recommended r-cran-*
RUN apt -y autoremove
RUN apt -y install software-properties-common
RUN apt update
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
RUN apt update
RUN apt -y install r-base r-base-core r-recommended r-base-dev

RUN apt-get -y install r-base-dev
#    && apt-get -y install r-cran-data.table\
#    && apt-get -y install r-bioc-biobase\
#    && apt-get -y install r-bioc-biocgenerics\
#    && apt-get -y install r-bioc-deseq2\
#    && apt-get -y install r-cran-randomforest\
#    && apt-get -y install r-cran-mvtnorm\
#    && apt-get -y install r-bioc-biocinstaller\
#    && apt-get -y install r-cran-ade4\
#    && apt-get -y install r-cran-minqa

RUN apt-get clean
RUN locale-gen en_US.UTF-8 en fi_FI.UTF-8

RUN mkdir /src


# INSTALL PYPROPHET
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py
RUN python2 get-pip.py
#RUN pip install numpy
RUN pip install pyprophet==0.24.1

# INSTALL TPP
RUN mkdir -p /opt/tpp/
RUN mkdir /opt/tpp-data
WORKDIR /src/
RUN svn checkout svn://svn.code.sf.net/p/sashimi/code/tags/release_5-2-0
RUN echo "INSTALL_DIR = /opt/tpp\nBASE_URL = /tpp\nTPP_DATADIR = /opt/tpp-data" > release_5-2-0/site.mk
WORKDIR /src/release_5-2-0
COPY tpp-5.2-fix.diff /root/
COPY comet_source_2019015.zip extern/
RUN cat /root/tpp-5.2-fix.diff |patch -p1
RUN make libgd
RUN make all
RUN make install

# INSTALL msproteomicstools.git
WORKDIR /src/
# RUN pip3 install msproteomicstools
WORKDIR /src/
RUN git clone https://github.com/msproteomicstools/msproteomicstools.git
WORKDIR /src/msproteomicstools
RUN git checkout v0.6.0
RUN python3 setup.py install

# Install Comet
RUN mkdir -p /opt/comet
WORKDIR /opt/comet
RUN wget https://sourceforge.net/projects/comet-ms/files/comet_2019015.zip
RUN unzip comet_2019015.zip
RUN ln -s comet.2019015.linux.exe comet-ms
RUN chmod ugo+x comet.2019015.linux.exe

# Install tandem
WORKDIR /opt
RUN wget ftp://ftp.thegpm.org/projects/tandem/source/tandem-linux-17-02-01-4.zip
RUN unzip tandem-linux-17-02-01-4.zip
RUN mv tandem-linux-17-02-01-4 tandem
RUN ln -s /opt/tandem/bin/static_link_ubuntu/tandem.exe /opt/tandem/tandem
RUN chmod ugo+x /opt/tandem/bin/static_link_ubuntu/tandem.exe

# INSTALL msgfplus
#RUN mkdir /opt/msgfplus
#WORKDIR /opt/msgfplus
#RUN wget https://github.com/MSGFPlus/msgfplus/releases/download/v2018.10.15/MSGFPlus_v20181015.zip
#RUN unzip MSGFPlus_v20181015.zip

# INSTALL Percolator
WORKDIR /opt
RUN wget https://github.com/percolator/percolator/releases/download/rel-3-01/ubuntu64_release.tar.gz
RUN tar xfv ubuntu64_release.tar.gz
RUN dpkg -i percolator-converters-v3-01-linux-amd64.deb percolator-v3-01-linux-amd64.deb

# INSTALL luciphor2
RUN mkdir /opt/luciphor2
WORKDIR /opt/luciphor2
RUN wget https://sourceforge.net/projects/luciphor2/files/luciphor2.jar

# INSTALL dia-umpire
RUN mkdir /opt/dia-umpire
WORKDIR /opt/dia-umpire
RUN wget https://github.com/guoci/DIA-Umpire/releases/download/v2.1.3/v2.1.3.zip
RUN unzip v2.1.3.zip
RUN ln -s v2.1.3/DIA_Umpire_SE.jar DIA_Umpire_SE.jar

## Fetch gladiator and install needed R-packages
RUN mkdir /opt/gladiator
COPY comet_settings_template.xml xtandem_settings_template.xml dia-pipeline.py install-R-packages.R iRTAssayLibrary.TraML iRT.txt diaumpire-params-template.txt swaths2stats.R /opt/gladiator/
RUN mkdir /.Rcache
RUN mkdir /opt/Rlibs/
RUN chmod u+x /opt/gladiator/install-R-packages.R
RUN /opt/gladiator/install-R-packages.R

WORKDIR /root

# UI
RUN pip3 install "pyramid==1.10.2" waitress
RUN pip3 install cookiecutter

WORKDIR /

# Install thermo raw library (disabled)
# COPY wineprefix64.tar.gz /
# RUN tar xfv wineprefix64.tar.gz
# RUN dpkg --add-architecture i386
# RUN wget -nc https://dl.winehq.org/wine-builds/winehq.key
# RUN sudo apt-key add winehq.key
# RUN apt-add-repository 'deb https://dl.winehq.org/wine-builds/ubuntu/ focal main'
# RUN apt -y install --install-recommends winehq-stable
# RUN ln -s /wineprefix64/ /root/.wine
# RUN apt -y install xvfb
# RUN apt -y install winetricks

ENV PATH=${PATH}:/ThermoRawFileParser
COPY install-ThermoRawFileParser.sh /root
RUN chmod u+x /root/install-ThermoRawFileParser.sh

WORKDIR /workdir
ENV PYTHONPATH=/opt/gladiator:/opt/gladiator/UI


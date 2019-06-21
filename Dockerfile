#FROM nvidia/cuda:latest
#FROM nvidia/cuda:8.0-cudnn7-devel-ubuntu16.04
FROM nvidia/cuda:10.0-cudnn7-runtime-ubuntu16.04

MAINTAINER chenaoki <chenaoki@gmail.com>

RUN echo "now building"

RUN apt-get update --fix-missing && apt-get -y install \
                   build-essential \
                   git vim curl wget zsh make ffmpeg wget bzip2 \
                   ca-certificates \
                   zlib1g-dev \
                   libssl-dev \
                   libbz2-dev \
                   libreadline-dev \
                   libsqlite3-dev \
                   libglib2.0-0 \
                   libxext6 \
                   libsm6 \
                   libxrender1 \
                   llvm \
                   libncurses5-dev \
                   libncursesw5-dev \
                   libpng-dev \
                   libgtk2.0-0 \
                   mercurial \
                   subversion \
                   python-qt4 


USER root
ENV HOME /root
ENV NOTEBOOK_HOME /notebooks
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

##########################
# Anaconda installation

RUN wget --quiet https://repo.anaconda.com/archive/Anaconda2-4.3.1-Linux-x86_64.sh -O ~/anaconda.sh && \
  /bin/bash ~/anaconda.sh -b -p /opt/conda && \
  rm ~/anaconda.sh && \
  ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
  echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
  echo "conda activate base" >> ~/.bashrc

RUN apt-get install -y curl grep sed dpkg && \
  TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
  curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
  dpkg -i tini.deb && \
  rm tini.deb && \
  apt-get clean

RUN conda update -y conda
RUN conda install -c numba numba=0.33.0
RUN conda install -c numba -y pyculib=1.0.1
RUN conda install --override-channels -c anaconda chainer 
#RUN conda install -c menpo opencv3


##########################

RUN pip install --upgrade pip
RUN pip install cupy-cuda100

##########################


WORKDIR $HOME 
RUN jupyter notebook --generate-config
RUN echo "c.NotebookApp.ip = '*'" >> $HOME/.jupyter/jupyter_notebook_config.py
RUN echo "c.NotebookApp.port = 8888" >> $HOME/.jupyter/jupyter_notebook_config.py
RUN echo "c.NotebookApp.open_browser = False" >> $HOME/.jupyter/jupyter_notebook_config.py
RUN echo "c.NotebookApp.notebook_dir = '$NOTEBOOK_HOME'" >> $HOME/.jupyter/jupyter_notebook_config.py


##########################


RUN echo "cloning dotfiles"
WORKDIR $HOME
RUN git clone https://github.com/chenaoki/dotfiles.git
WORKDIR $HOME/dotfiles
RUN python install.py

##########################


RUN mkdir -p $NOTEBOOK_HOME 
CMD ["sh", "-c", "jupyter notebook > $NOTEBOOK_HOME/log.txt 2>&1"]

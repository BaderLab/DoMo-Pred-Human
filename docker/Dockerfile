FROM python:2.7.18

#install required packages 

#install Cython
RUN pip install Cython==0.22

# instal setuptools
RUN pip install setuptools==17.1.1

#install numpy
RUN pip install numpy==1.9.2

#install scipy
RUN pip install scipy==0.16.0

#copy the domo_pred code and files into the image
COPY DoMo_Pred /DoMo_Pred

#install additional packages that weren't listed as requirements but 
# are still needed for compilation
RUN pip install networkx==1.11
RUN pip install sklearn==0.0
RUN pip install dill

#compile the PWMsearch
WORKDIR /DoMo_Pred/Peptide/PWMsearch/source
RUN g++ -fPIC -c PWMSearch.cpp -I/usr/local/include/python2.7 -lpython2.7 -o PWMSearch.o
RUN g++ -shared PWMSearch.o -o PWMSearch.so

#install the NW align library
WORKDIR /DoMo_Pred/Peptide/Structure/source/NW/nwalign/nwalign-0.3.1
RUN python2.7 setup.py install

#install vim
RUN apt-get update
RUN apt-get install -y vim

WORKDIR /DoMo_Pred

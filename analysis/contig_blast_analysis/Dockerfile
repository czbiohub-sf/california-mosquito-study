FROM continuumio/anaconda3
MAINTAINER lucy.li@czbiohub.org

USER root

# Add user "main" because that's what is expected by this image
RUN useradd -ms /bin/bash main

RUN pip install biopython
RUN pip install ete3
RUN pip install pandas
RUN pip install boto3

ADD lca_analysis.py /usr/local/bin/lca_analysis.py
ADD lca_functions.py /usr/local/bin/lca_functions.py

RUN chmod +x /usr/local/bin/lca_analysis.py
RUN chmod +x /usr/local/bin/lca_functions.py

WORKDIR /home

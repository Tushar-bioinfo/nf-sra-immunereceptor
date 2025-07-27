FROM continuumio/miniconda3

COPY environment.yml /tmp/environment.yml

RUN conda env create -f /tmp/environment.yml

SHELL ["conda", "run", "--no-capture-output", "-n", "sra", "/bin/bash", "-c"]
 #Set working directory 
WORKDIR /workspace
# Set default shell 
ENTRYPOINT ["/bin/bash"]

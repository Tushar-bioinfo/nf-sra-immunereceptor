
FROM continuumio/miniconda3

COPY environment.yml /tmp/environment.yml

RUN conda env create -f /tmp/environment.yml

SHELL ["conda", "run", "-n", "immune-receptor-pipeline", "/bin/bash", "-c"]

# Set entry point to bash
ENTRYPOINT ["/bin/bash"]

# syntax=docker/dockerfile:1

# A Dockerfile that can be used without cloning the project manually. It uses
# malb/latticeestimator at HEAD. This is appropriate for users who want to use
# malb/latticeestimator but do not want to set up a local development
# environment.
FROM sagemath/sagemath:latest

user root
RUN apt update && apt upgrade -y
RUN apt install -y git

RUN git clone https://github.com/malb/lattice-estimator.git /latticeestimator
WORKDIR /lattice-estimator

RUN sage -pip install -r /latticeestimator/requirements.txt
CMD ["sage"]
